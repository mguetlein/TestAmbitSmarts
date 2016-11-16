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
	// https://sourceforge.net/p/ambit/bugs/109/
	public void smirksIncludingDegreeOfConnectionsNotWorking_4297_c0857()
	{
		String smirks_4297_old = "[#8:4]-,=[#6;D3R0:1](-[#6:2])!@-[#6:3]-[#6:8](=[O:10])-[#6:7](-[#8-:6])=[O:9]>>[#6:3]-[#6:8](=[O:10])-[#6:7](-[#8-:6])=[O:9].[#6:2]-[#6:1](-[#8-])=[O:4]";
		String smirks_4297_new = "[#8:4]([H,!H])-,=[#6;D3,D4H1;R0:1]([H,!H])(-[#6:2])!@-[#6:3]-[#6:8](=[O:10])-[#6:7](-[#8-:6])=[O:9]>>[#6:3]-[#6:8](=[O:10])-[#6:7](-[#8-:6])=[O:9].[#6:2]-[#6:1](-[#8-])=[O:4]";

		String expectedSmiles1;
		String expectedSmiles2;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("CC(=O)C(=O)[O-]");
			expectedSmiles1 = SmilesGenerator.absolute().create(mol);
			mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("C(C(=O)[O-])Cl");
			expectedSmiles2 = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		System.out.println(expectedSmiles1);
		System.out.println(expectedSmiles2);

		for (String smirks : new String[] { smirks_4297_old, smirks_4297_new })
		{
			String smi_c1248 = "OC(CCl)CC(=O)C([O-])=O";
			List<String> s = applySmirks(smirks, smi_c1248);

			if (smirks.equals(smirks_4297_old))
			{
				Assert.assertEquals(s.size(), 0);
			}
			else
			{
				Assert.assertTrue(s.get(0).contains(expectedSmiles1));
				Assert.assertTrue(s.get(0).contains(expectedSmiles2));
			}
		}
	}

	/**
	 * the rings with the =O are actually not aromatic
	 */
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
		// change to false, as it should not match
		Assert.assertFalse(s.contains(expectedSmiles));
	}

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

	public static void main(String[] args)
	{
		TestAmbitWithdrawn t = new TestAmbitWithdrawn();
		t.smirksIncludingDegreeOfConnectionsNotWorking_4297_c0857();
	}
}
