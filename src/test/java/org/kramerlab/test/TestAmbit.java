package org.kramerlab.test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.core.data.MoleculeTools;
import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.SmartsConst;

@RunWith(JUnit4.class)
public class TestAmbit
{
	@Test
	// https://sourceforge.net/p/ambit/bugs/108/
	public void aromBondNotAddedInSmirksWithSmartsBonds_rule2690_c0138()
	{
		String smirks = "[#6:9]=,:[#6:10]@-[C:6]([H])([#8:2][H:11])[C:5]([H])(@-[#6:7]=,:[#6:8])[#8:1][H:12]>>[H:12][#8:1]-[c:5](:[c:7]:[c:8]):[c:6](-[#8:2][H:11]):[c:10]:[c:9]";
		String smi = "O[C@@H]1C=CC=C(Cl)[C@@H]1O";
		String expectedSmiles;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("C1=CC(=C(C(=C1)Cl)O)O");
			expectedSmiles = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue(s.contains(expectedSmiles));
	}

	@Test
	// https://sourceforge.net/p/ambit/bugs/106/
	public void smirksNotMatchingAllRingsInPolycyclicAromaticCompound_rule4282_c0857()
	{
		String smirks = "[H][c:1]1[c:7][c:10][c:9][#6,#7;a:8][c:2]1[H]>>[#8]-[c:1]1[c:7][c:10][c:9][#6,#7;a:8][c:2]1-[#8]";
		String smi = "O=c1ccc2ccc3ccc(=O)c4ccc1c2c34";
		String expectedSmiles;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("C1=CC2=C3C(=CC=C4C(=O)C=CC1=C43)C(=O)C(=C2O)O");
			expectedSmiles = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue(s.contains(expectedSmiles));
	}

	@Test
	// https://sourceforge.net/p/ambit/bugs/107/
	public void chiralityInformationNotAdded_rule2978_c0105()
	{
		String smirks = "[c:1]1[c:6]([H])[c:5]([H])[c:4][c:3][c:2]1>>[#8]([H])-[#6@H:5]-1-[#6:4]=[#6:3]-[#6:2]=[#6:1]-[#6@H:6]-1-[#8]([H])";
		String smi = "Clc1ccccc1";
		String expectedSmiles;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("C1=C[C@H]([C@H](C(=C1)Cl)O)O");
			expectedSmiles = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue(s.contains(expectedSmiles));
	}

	public void ringSplit2_rule4287()
	{
		String smirks = "[#8-:10]-[#6:9](=[O:11])-[#6:1]-1=[CH1]-[#6;H1:5]=[#6;H1:4]-[#6;H1:3]=[#6;H1:2]-1>>O-[#6:5](=O)-[#6;H2:4]\\[#6;H1:3]=[#6;H1:2]/[#6;H2:1]-[#6:9](-[#8-:10])=[O:11]";
		String smi = "O=C([O-])C1=CC=CC=C1";
		List<String> products = applySmirks(smirks, smi);
		Map<String, String> diff = diffProductToReference(products,
				new String[] { "OC(=O)C\\C=C/CC([O-])=O" });
		Assert.assertTrue("" + diff, null == diff);
	}

	@Test
	public void ringSplit3_rule3743_u50144_u137948()
	{
		String smirks = "[#6:3]-[#8:15]-[c;R:11]1[c:12](-[#8:2]([H]))[c:7](-[#8:1]([H]))[c;R:8]([#1,#6,#9,#17,#35,#53;A:4])[c;R:9](-[!#16:6])[c;R:10]1-[!#8!#16:5]>>[!#16:6]\\[#6:9](=[#6:10](/[!#8!#16:5])-[#6:11](-[#8:15])=O)-[#6:8]([#1,#6,#9,#17,#35,#53;A:4])-[#6:7](=[O:1])-[#6:12](-[#8-:2])=O.[#6:3]-[#8]";
		String smi1 = "COC1=C2C=CC=CC2=CC(=C1O)O";
		String smi2 = "C1=CC2=CC(=C(C3=C2C(=C1)CC(=O)O3)O)O";
		Assert.assertFalse("Results should not be empty",
				applySmirks(smirks, smi1).isEmpty() || applySmirks(smirks, smi2).isEmpty());
	}

	@Test
	public void ringSplit4_rule4298_u140722()
	{
		String smirks = "[#6:11]@-[c:2]1[c;R1:4][c;R1:5][c:6](-[#8:8]([H]))[c:7](-[#8:1]([H]))[c:3]1@-[#6:12]>>[#6:12]-[#6:3](=O)-[#6:2](\\[#6:11])=[#6:4]/[#6:5]=[#6:6](/[#8:8]([H]))-[#6:7](-[O-])=[O:1]";

		String onceAromatic = "Oc1ccc2CCCCc2c1O";
		List<String> products = applySmirks(smirks, onceAromatic);
		Map<String, String> diff = diffProductToReference(products,
				new String[] { "O\\C(=C\\C=C1\\CCCCC1=O)C([O-])=O" });
		Assert.assertTrue("" + diff, null == diff);

		String twiceAromatic = "Oc1ccc2ccccc2c1O";
		products = applySmirks(smirks, twiceAromatic);
		diff = diffProductToReference(products,
				new String[] { "O\\C(=C\\C=C1\\C=CC=CC1=O)C([O-])=O" });
		Assert.assertTrue("" + diff, null == diff);

		String smiles_u140722 = "OCc1cc(O)c(O)c2cccc(C([O-])=O)c12";
		products = applySmirks(smirks, smiles_u140722);
		diff = diffProductToReference(products,
				new String[] { "C1=CC(=O)\\C(=C(/C=C(\\C(=O)[O-])/O)\\CO)\\C(=C1)C(=O)[O-]" });
		Assert.assertTrue("" + diff, null == diff);
	}

	@Test
	public void ringSplit5_rule4224_missing()
	{
		String smirks = "[#8:2]([H])-[c:10]1[#6,#7;a:9][c:8][c:7][c:6](-[*,#1:11])[c:5]1-[#8:1]([H])>>[#8-:2]-[#6:10](=O)\\[#6,#7:9]=[#6:8]/[#6:7]=[#6:6](/[*,#1:11])-[#6:5](-[#8-:1])=O";

		String u97336 = "Oc1ncc2ccccc2c1O";
		List<String> products = applySmirks(smirks, u97336);
		Map<String, String> diff = diffProductToReference(products,
				new String[] { "C1=CC(=C(C=C1)C(=O)[O-])/C=N\\C(=O)[O-]" });
		Assert.assertTrue("" + diff, null == diff);

		String u8723 = "CC(C(C)=O)c1cc2c(Cl)nc(O)c(O)c2[nH]1";
		products = applySmirks(smirks, u8723);
		diff = diffProductToReference(products,
				new String[] { "CC(C(C)=O)c1cc(\\C(Cl)=N/C([O-])=O)c([nH]1)C([O-])=O" });
		Assert.assertTrue("" + diff, null == diff);
	}

	@Test
	public void ringSplit5_rule4224_toomany()
	{
		String smirks = "[#8:2]([H])-[c:10]1[#6,#7;a:9][c:8][c:7][c:6](-[*,#1:11])[c:5]1-[#8:1]([H])>>[#8-:2]-[#6:10](=O)\\[#6,#7:9]=[#6:8]/[#6:7]=[#6:6](/[*,#1:11])-[#6:5](-[#8-:1])=O";

		String smiles = "Nc1ccc(O)c(O)c1";
		List<String> products = applySmirks(smirks, smiles);
		Map<String, String> diff = diffProductToReference(products,
				new String[] { "N\\C(\\C=C/C([O-])=O)=C\\C([O-])=O" });
		Assert.assertTrue("" + diff, null == diff);

	}

	@Test
	public void ringSplit5_rule4224_valenceError()
	{
		String smirks = "[#8:2]([H])-[c:10]1[#6,#7;a:9][c:8][c:7][c:6](-[*,#1:11])[c:5]1-[#8:1]([H])>>[#8-:2]-[#6:10](=O)\\[#6,#7:9]=[#6:8]/[#6:7]=[#6:6](/[*,#1:11])-[#6:5](-[#8-:1])=O";

		String c0707 = "OCc1ccc2ccc(O)c(O)c2c1";
		List<String> products = applySmirks(smirks, c0707);
		Map<String, String> diff = diffProductToReference(products,
				new String[] { "C1=C(C=C(C(=C1)/C=C\\C(=O)[O-])C(=O)[O-])CO" });
		Assert.assertTrue("" + diff, null == diff);

	}

	public static Map<String, String> diffProductToReference(List<String> products,
			String[] reference)
	{
		Set<String> p = new HashSet<String>();
		p.addAll(products);
		Set<String> r = new HashSet<String>();
		r.addAll(Arrays.asList(reference));
		diffSmiles(p, r);
		if (p.isEmpty() && r.isEmpty())
			return null;
		Map<String, String> diff = new HashMap<String, String>();
		for (String notInReference : p)
		{
			diff.put(notInReference, "notInReference");
		}
		for (String notInProducts : r)
		{
			diff.put(notInProducts, "notInProducts");
		}
		return diff;
	}

	public static List<String> applySmirks(String smrk, String smi)
	{
		try
		{
			SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
			smrkMan.setFlagSSMode(SmartsConst.SSM_MODE.SSM_NON_IDENTICAL_FIRST);
			smrkMan.setFlagProcessResultStructures(true);
			smrkMan.setFlagClearHybridizationBeforeResultProcess(true);
			smrkMan.setFlagClearImplicitHAtomsBeforeResultProcess(true);
			smrkMan.setFlagClearAromaticityBeforeResultProcess(true);
			smrkMan.setFlagAddImplicitHAtomsOnResultProcess(true);
			smrkMan.setFlagConvertAddedImplicitHToExplicitOnResultProcess(false);
			smrkMan.setFlagConvertExplicitHToImplicitOnResultProcess(true);
			smrkMan.getSmartsParser().mSupportDoubleBondAromaticityNotSpecified = false;
			smrkMan.setFlagApplyStereoTransformation(true);

			SMIRKSReaction reaction = smrkMan.parse(smrk);
			if (!smrkMan.getErrors().equals(""))
				throw new RuntimeException("Invalid SMIRKS: " + smrkMan.getErrors());

			IAtomContainer target = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles(smi);
			for (IAtom atom : target.atoms())
				if (atom.getFlag(CDKConstants.ISAROMATIC))
					atom.setFlag(CDKConstants.ISAROMATIC, false);
			for (IBond bond : target.bonds())
				if (bond.getFlag(CDKConstants.ISAROMATIC))
					bond.setFlag(CDKConstants.ISAROMATIC, false);
			// do not add Hs: https://sourceforge.net/p/cdk/mailman/message/34608714/
			//			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);
			//			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
			//			adder.addImplicitHydrogens(target);
			//CDK bug regarding stereo in 1.5.13
			//AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);
			MoleculeTools.convertImplicitToExplicitHydrogens(target);

			// for our project, we want to use daylight aromaticity
			// CDKHueckelAromaticityDetector.detectAromaticity(target);
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(),
					Cycles.or(Cycles.all(), Cycles.edgeShort()));
			aromaticity.apply(target);

			IAtomContainerSet resSet2 = smrkMan.applyTransformationWithSingleCopyForEachPos(target,
					null, reaction, SmartsConst.SSM_MODE.SSM_ALL);
			List<String> result = new ArrayList<String>();
			if (resSet2 != null)
				for (int i = 0; i < resSet2.getAtomContainerCount(); i++)
				{
					IAtomContainer mol = resSet2.getAtomContainer(i);
					AtomContainerManipulator.suppressHydrogens(mol);
					String smiles = SmilesGenerator.absolute().create(mol);
					result.add(smiles);
				}
			return result;
		}
		catch (Exception e)
		{
			System.err.println("SMIRKS " + smrk);
			System.err.println("SMILES " + smi);
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * compares smiles string
	 * @param smiA
	 * @param smiB
	 * @return true if smiA and smiB are equal
	 */
	public static boolean sameSame(String smiA, String smiB)
	{
		/*
		 * TODO: compare standardized smiles
		 */
		return smiA.equals(smiB);
	}

	/**
	 * removes intersection from both sets, based on the {@link #sameSame(String, String)} comparing method
	 * @param smilesA
	 * @param smilesB
	 */
	public static void diffSmiles(Set<String> smilesA, Set<String> smilesB)
	{
		for (String smiA : smilesA)
		{
			for (String smiB : smilesB)
			{
				if (sameSame(smiA, smiB))
				{
					smilesA.remove(smiA);
					smilesB.remove(smiB);
					diffSmiles(smilesA, smilesB);
					return;
				}
			}
		}
	}

	public static void main(String[] args) throws CDKException
	{
		TestAmbit t = new TestAmbit();
		t.aromBondNotAddedInSmirksWithSmartsBonds_rule2690_c0138();
		//		t.productToSmilesError_rule3707_u143203();
	}
}
