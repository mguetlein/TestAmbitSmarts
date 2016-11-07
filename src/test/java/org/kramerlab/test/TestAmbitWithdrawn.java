package org.kramerlab.test;

import java.util.List;

import org.junit.Assert;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

@RunWith(JUnit4.class)
public class TestAmbitWithdrawn
{
	// https://sourceforge.net/p/ambit/bugs/108/
	public void aromBondNotAddedInSmirksWithSmartsBonds_rule2690_c0138()
	{
		/**
		 * ambit cannot introduce aromatic form out of not-arom, instead =- bonds have to be introduced 
		 * https://github.com/enviPath/enviPath/issues/101
		 */
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
