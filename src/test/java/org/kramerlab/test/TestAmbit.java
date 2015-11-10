package org.kramerlab.test;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.SmartsConst;

@RunWith(JUnit4.class)
public class TestAmbit
{
	@Test
	public void stereoChemIsLost_rule4212_u136348()
	{
		String smirks = "[H:6][C:1]([#6:4])([#16;H1v2])[#1,#6:5]>>[H:6][C:1]([H])([#6:4])[#1,#6:5]";
		String smi = "CN\\C(NCCS)=C\\[N+]([O-])=O";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should still contain stereocheminfo: " + s.get(0), s.get(0).contains("\\"));
	}

	@Test
	public void ringSplitArom_rule4224_ar13()
	{
		String smirks = "[#8;H1:2]-[#6:10]1:[#6:5](-[#8;H1:1]):[#6:6](-[*,#1:11]):[#6:7]:[#6:8]:[#6,#7:9]:1>>[#8&-:2]-[#6:10](=O)-[#6,#7:9]=[#6:8]-[#6:7]=[#6:6](-[*,#1:11])-[#6:5](-[#8&-:1])=O";
		String smi = "O-C(=C-C=C-C=O)C([O-])=O";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	@Test
	public void ringSplitKekulized_rule4224_ar13()
	{
		String smirks = "[#8;H1:2]-[#6:10]1=[#6:5](-[#8;H1:1])-[#6:6](-[*,#1:11])=[#6:7]-[#6:8]=[#6,#7:9]-1>>[#8&-:2]-[#6:10](=O)-[#6,#7:9]=[#6:8]-[#6:7]=[#6:6](-[*,#1:11])-[#6:5](-[#8&-:1])=O";
		String smi = "O-C(=C-C=C-C=O)C([O-])=O";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	@Test
	public void anotherRingSplit_rule4287()
	{
		String smirks = "[#8-:10]-[#6:9](=[O:11])-[#6:1]-1=[CH1]-[#6;H1:5]=[#6;H1:4]-[#6;H1:3]=[#6;H1:2]-1>>O-[#6:5](=O)-[#6;H2:4]\\[#6;H1:3]=[#6;H1:2]/[#6;H2:1]-[#6:9](-[#8-:10])=[O:11]";
		String smi = "O=C([O-])C1=CC=CC=C1";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	@Test
	public void complexCarbonStructure_rule3803_c0630()
	{
		// was working in cdk1.4
		String smirks = "[H][C;R1:3]([H:2])([#6;X4:5])[#6;X4:6]>>[H:2][C;R1:3]([#6;X4:5])([#6;X4:6])O";
		String smi = "CC1(C)C2CC1C(CO)=CC2";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertFalse("Results should not be empty", s.isEmpty());
	}

	@Test
	public void ringSplitHeteroAtom_rule4294_u78282()
	{
		String smirks = "[O:8]=[c:2]1[c:3][c:4][c:5][c:6][o:1]1>>[#8:1]([H])\\[#6:6]=[#6:5]/[#6:4]=[#6:3]\\[#6:2](-[#8-])=[O:8]";
		String smi = "C1=CC=C2C(=C1)C=CC(=O)O2";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	@Test
	public void missingProducts_rule2793_u26103()
	{
		String smirks = "[#8:8]([H])-[c:2]1[c:1](-[#8:7]([H]))[c;R]([c;R:5](-[!#8,#1:11])[c;R:4](-[!#8,#1:10])[c;R:3]1-[!#8,#1:9])S([#8])(=O)=O>>[!#8,#1:11]\\\\\\[#6:5]=[#6:4](///[!#8,#1:10])-[#6:3](-[!#8,#1:9])-[#6:2](=[O:8])-[#6:1](-[#8-])=[O:7]";
		String smi = "C1=C(C=CC(=C1)C2=NC3=C(C(=C(C=C3N2)S(=O)(=O)[O-])O)O)C4=NC5=C(C=C(C=C5S(=O)(=O)[O-])S(=O)(=O)[O-])N4";
		// working smiles example: C1=CC(=C(C(=C1N)S(=O)(=O)[O-])O)O
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	@Test
	public void retainUnfilledValence_rule1196_u133400()
	{
		// the unfilled valence of nitrogen (-> [N]) should remaine unchanged
		String smirks = "[H][#6:1](-[#6:5])=[O:4]>>[#6:5]-[#6:1](-[#8-])=[O:4]";
		String smi = "C[N]C(=O)C(=O)C=O";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should still contain [N]: " + s.get(0), s.get(0).contains("[N]"));

		// however, at the same time, Hs should be added to newly created atoms  
		smirks = "[H:5][C:1]([#6:6])([#1,#9,#17,#35,#53:4])[#9,#17,#35,#53]>>[H:5][C:1]([#6:6])([#8])[#1,#9,#17,#35,#53:4]";
		smi = "C(CN(CCCl)CC(C(=O)O)N)Cl";
		s = applySmirks(smirks, smi);
		Assert.assertFalse("Results should NOT contain [O]: " + s.get(0), s.get(0).contains("[O]"));
	}

	// fixed by changing the smirks

	@Test
	public void oneHToMany_rule3769_c0107()
	{
		//String smirks = "[#8H1:7]-[#6:6](-[#6:8](-[#8-:9])=[O:10])=[#6:5](-[#1,#6,#17:11])-[#6:1]=[#6:2]-[#6;R0:3](-[#1,#6,#16:13])=[O:12]>>[#8-:9]-[#6:8](=[O:10])-[#6:6](=[O:7])-[#6:5](-[#1,#6,#17:11])-[#6:1]=[#6:2].[O-]-[#6:3](-[#1,#6,#16:13])=[O:12]";
		String smirks = "[H][#8:7]-[#6:6](-[#6:8](-[#8-:9])=[O:10])=[#6:5](-[#1,#6,#17:11])-[#6:1]=[#6:2]-[#6;R0:3](-[#1,#6,#16:13])=[O:12]>>[#8-:9]-[#6:8](=[O:10])-[#6:6](=[O:7])-[#6:5](-[#1,#6,#17:11])-[#6:1]=[#6:2].[O-]-[#6:3](-[#1,#6,#16:13])=[O:12]";
		String smi = "O-C(=C-C=C-C=O)C([O-])=O";
		List<String> s = applySmirks(smirks, smi);
		for (String p : s)
			Assert.assertFalse("should not contain [OH]= : '" + p + "'", p.contains("[OH]="));
	}

	// cis/trans smirks do not throw error anymore

	@Test
	public void smirksWithCisTrans_rule3908()
	{
		String smirks = "[#8;H1:2]-[c:12]1[c:7](-[#8;H1:1])[c;R1:8]([#1,#6,#9,#17,#35,#53;A:3])[c;R1:9](-[!#16:6])[c;R1:10](-[!#8!#16:5])[c;R1:11]1-[#1,#6,#7:4]>>[#8-:2]-[#6:12](=O)-[#6:7](=[O:1])-[#6:8]([#1,#6,#9,#17,#35,#53;A:3])-[#6:9](\\[!#16:6])=[#6:10]\\[!#8!#16:5].[O-]-[#6:11](-[#1,#6,#7:4])=O";
		String smi = "C";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertNotNull("SMIRKS parsing should not fail", s);
	}

	// working, problem was non-kekulized input
	@Test
	public void yetAnotherRingSplit_rule4185_u132707()
	{
		String smirks = "[#8:7]([H])-[c:2]1[c:6]([H])[c:5]([H])[c:4](-[#8:8]([H]))[n:1][c:3]([H])1>>[#8-:7]-[#6:2](=O)[#6:6]=[#6:5]/[#6:4](=[O:8])-[#7:1]-[#6:3]=O";
		String smi = "C1=CC(=NC=C1O)O";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should not contain C(=C=C: " + s.get(0), s.get(0).contains("C(=C=C"));
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

			SMIRKSReaction reaction = smrkMan.parse(smrk);
			if (!smrkMan.getErrors().equals(""))
				throw new RuntimeException("Invalid SMIRKS: " + smrkMan.getErrors());

			IAtomContainer target = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smi);
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
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);

			// for our project, we want to use daylight aromaticity
			// CDKHueckelAromaticityDetector.detectAromaticity(target);
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(),
					Cycles.edgeShort()));
			aromaticity.apply(target);

			IAtomContainerSet resSet2 = smrkMan.applyTransformationWithSingleCopyForEachPos(target, null, reaction,
					SmartsConst.SSM_MODE.SSM_ALL);
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

	public static void main(String[] args)
	{
		TestAmbit t = new TestAmbit();
		//		t.missingProducts_rule2793_u26103();
		t.retainUnfilledValence_rule1196_u133400();
	}
}
