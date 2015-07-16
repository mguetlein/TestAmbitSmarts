package org.kramerlab.test;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.core.helper.CDKHueckelAromaticityDetector;
import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.SmartsConst;
import ambit2.smarts.SmartsHelper;

@RunWith(JUnit4.class)
public class TestAmbit
{
	@Test
	public void oneHToMany_rule3769_c0107()
	{
		String smirks = "[#8H1:7]-[#6:6](-[#6:8](-[#8-:9])=[O:10])=[#6:5](-[#1,#6,#17:11])-[#6:1]=[#6:2]-[#6;R0:3](-[#1,#6,#16:13])=[O:12]>>[#8-:9]-[#6:8](=[O:10])-[#6:6](=[O:7])-[#6:5](-[#1,#6,#17:11])-[#6:1]=[#6:2].[O-]-[#6:3](-[#1,#6,#16:13])=[O:12]";
		String smi = "O-C(=C-C=C-C=O)C([O-])=O";
		List<String> s = applySmirks(smirks, smi);
		for (String p : s)
			Assert.assertFalse("should not contain [OH]= : '" + p + "'", p.contains("[OH]="));
	}

	@Test
	void ringSplitArom_rule4224_ar13()
	{
		String smirks = "[#8;H:2]-[#6:10]1:[#6:5](-[#8;H:1]):[#6:6](-[*,#1:11]):[#6:7]:[#6:8]:[#6,#7:9]:1>>[#8&-:2]-[#6:10](=O)-[#6,#7:9]=[#6:8]-[#6:7]=[#6:6](-[*,#1:11])-[#6:5](-[#8&-:1])=O";
		String smi = "O-C(=C-C=C-C=O)C([O-])=O";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	@Test
	void ringSplitKekulized_rule4224_ar13()
	{
		String smirks = "[#8;H:2]-[#6:10]1=[#6:5](-[#8;H:1])-[#6:6](-[*,#1:11])=[#6:7]-[#6:8]=[#6,#7:9]-1>>[#8&-:2]-[#6:10](=O)-[#6,#7:9]=[#6:8]-[#6:7]=[#6:6](-[*,#1:11])-[#6:5](-[#8&-:1])=O";
		String smi = "O-C(=C-C=C-C=O)C([O-])=O";
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

	// cis/trans smirks do not throw error anymore

	@Test
	public void smirksWithCisTrans_rule3908()
	{
		String smirks = "[#8;H1:2]-[c:12]1[c:7](-[#8;H1:1])[c;R1:8]([#1,#6,#9,#17,#35,#53;A:3])[c;R1:9](-[!#16:6])[c;R1:10](-[!#8!#16:5])[c;R1:11]1-[#1,#6,#7:4]>>[#8-:2]-[#6:12](=O)-[#6:7](=[O:1])-[#6:8]([#1,#6,#9,#17,#35,#53;A:3])-[#6:9](\\[!#16:6])=[#6:10]\\[!#8!#16:5].[O-]-[#6:11](-[#1,#6,#7:4])=O";
		String smi = "C";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertNotNull("SMIRKS parsing should not fail", s);
	}

	public static List<String> applySmirks(String smrk, String smi)
	{
		try
		{
			SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
			smrkMan.setFlagSSMode(SmartsConst.SSM_NON_IDENTICAL_FIRST);
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

			IAtomContainer target = SmartsHelper.getMoleculeFromSmiles(smi);
			for (IAtom atom : target.atoms())
				if (atom.getFlag(CDKConstants.ISAROMATIC))
					atom.setFlag(CDKConstants.ISAROMATIC, false);
			for (IBond bond : target.bonds())
				if (bond.getFlag(CDKConstants.ISAROMATIC))
					bond.setFlag(CDKConstants.ISAROMATIC, false);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
			adder.addImplicitHydrogens(target);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);
			CDKHueckelAromaticityDetector.detectAromaticity(target);

			IAtomContainerSet resSet2 = smrkMan.applyTransformationWithSingleCopyForEachPos(target, null, reaction,
					SmartsConst.SSM_ALL);
			List<String> result = new ArrayList<String>();
			if (resSet2 != null)
				for (int i = 0; i < resSet2.getAtomContainerCount(); i++)
					result.add(SmartsHelper.moleculeToSMILES(resSet2.getAtomContainer(i), true));
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
		t.ringSplitKekulized_rule4224_ar13();
	}
}
