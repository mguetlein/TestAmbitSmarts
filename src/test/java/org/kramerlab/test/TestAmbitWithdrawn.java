package org.kramerlab.test;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

@RunWith(JUnit4.class)
public class TestAmbitWithdrawn
{
	@Test
	public void complexCarbonStructure_rule3803_c0630()
	{
		/* was not empty in cdk1.4
		 * however: this is correct.
		 * R1 means exactly one ring and here we have a bridged molecule with two rings
		 * although somewhat hard to recognize.
		*/
		String smirks = "[H][C;R1:3]([H:2])([#6;X4:5])[#6;X4:6]>>[H:2][C;R1:3]([#6;X4:5])([#6;X4:6])O";
		String smi = "CC1(C)C2CC1C(CO)=CC2";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should be empty", s.isEmpty());
	}


	public static List<String> applySmirks(String smrk, String smi)
	{
		return TestAmbit.applySmirks(smrk, smi);
	}
}
