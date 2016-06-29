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
	public void ringSplit2_rule4287()
	{
		String smirks = "[#8-:10]-[#6:9](=[O:11])-[#6:1]-1=[CH1]-[#6;H1:5]=[#6;H1:4]-[#6;H1:3]=[#6;H1:2]-1>>O-[#6:5](=O)-[#6;H2:4]\\[#6;H1:3]=[#6;H1:2]/[#6;H2:1]-[#6:9](-[#8-:10])=[O:11]";
		String smi = "O=C([O-])C1=CC=CC=C1";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
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

	/*
	 * the problem is in the diverging notion of '@-'
	 * this is interpreted by Ambit as a bond in a ring, single or aromatic
	 * Chemaxon sees it as a bond in a ring, strictly single
	 * (cf. the notation for bonds without the '@' sign as in 2793:
	 * for Ambit '-' is strictly single, for Chemaxon it might also be aromatic)
	 */
	@Test
	public void aromaticRingWithOxygen_rule4150_u56188()
	{
		String smirks = "[#6:4]@-!:[#6;!$(C1(=O)C=CC(=O)C=C1)!$(C(=O)CC=O):1](@-!:[#6:2])=[O:5]>>[#6:2]@-[#8]@-[#6:1](@-[#6:4])=[O:5]";
		// simple true target
		String smi = "O=C1CCOC=C1";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertFalse("Results should not be empty", s.isEmpty());
		// simple false target: no ring
		smi = "O=C(C)C";
		s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should be empty", s.isEmpty());
		// simple false target: aromatic
		smi = "O=c1ccocc1";
		s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should be empty", s.isEmpty());
		// kekulized
		smi = "O=C1C=COC=C1";
		s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should be empty", s.isEmpty());
		// u56188 reduced to one target group only
		smi = "COc1cc(O)c2c(c1)oc(cc2=O)-c1ccc(O)cc1";
		s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should be empty", s.isEmpty());
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

	public static void main(String[] args) throws CDKException
	{
		TestAmbit t = new TestAmbit();
		t.ringSplit3_rule3743_u50144_u137948();
		//		t.productToSmilesError_rule3707_u143203();
	}
}
