package org.kramerlab.test;

import java.util.List;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

@RunWith(JUnit4.class)
public class TestAmbitFixed
{
	// fixed with stereo bug fix 14.06.2016
	@Test
	public void stereoChemNotInserted1_rule4230_u145861()
	{
		String smirks = "[#8-:11]-[#6:9](=[O:10])-[#6:1]([H])([H])-[c:2]1[c:3]([H])[c:4]([H])[c:5]([H])[c:6]([H])[c:7]([H])1>>[#8-:11]-[#6:9](=[O:10])\\[#6:1]([H])=[#6:2]-1\\[#6:7]([H])([H])-[#6:6]([H])=[#6:5]([H])-[#6:4]([H])=[#6:3]([H])-[#8]-1";
		String smi = "C1=CC=C(C=C1)CC(=O)[O-]";
		String expectedSmiles;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("[O-]C(=O)\\C=C1\\CC=CC=CO1");
			expectedSmiles = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		List<String> s = applySmirks(smirks, smi);
		Assert.assertEquals(expectedSmiles, s.get(0));
	}

	@Test
	public void stereoChemNotInserted2_rule2844_u114856()
	{
		String smirks = "[H][C:2]([#6:5]([H])([H])([H]))([#1,#6:4])!@-[#6:1]([H])([H])-[#6:3](-[#8-:8])=[O:6]>>[#6:5]([H])([H])([H])\\[#6:2](-[#1,#6:4])!@=[#6:1]\\[#6:3](-[#8-:8])=[O:6]";
		String smi = "CCC(C)CC(=O)[O-]";
		String expectedSmiles;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("CC\\C(C)=C/C([O-])=O");
			expectedSmiles = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		List<String> s = applySmirks(smirks, smi);
		Assert.assertEquals(expectedSmiles, s.get(0));
	}

	// fixed by rewriting smirks
	@Test
	public void makeRingAromatic_rule3667_u9685()
	{
		//String smirks = "[#8:7]([H])-[#6:1]([H])-1-[#6:2]=[#6:3]-[#6:4]=[#6:5]-[#6:6]([H])-1-[#8:8]([H])>>[#8:7]([H])-[c:1]1[c:2][c:3][c:4][c:5][c:6]1-[#8:8]([H])";
		String smirks = "[#8:7]([H])-[#6:1]([H])-1-[#6:2]=[#6:3]-[#6:4]=[#6:5]-[#6:6]([H])-1-[#8:8]([H])>>[#8:7]([H])-[#6:1]=1-[#6:2]=[#6:3]-[#6:4]=[#6:5]-[#6:6]=1-[#8:8]([H])";
		String smi = "C1=C[C@@H]([C@@H](C(=C1)C2=CC=C(C=C2)Cl)O)O";
		String expectedSmiles;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("Oc1cccc(c1O)-c1ccc(Cl)cc1");
			expectedSmiles = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		List<String> s = applySmirks(smirks, smi);
		Assert.assertEquals("Result be aromatic", expectedSmiles, s.get(0));
	}

	// fixed by using daylight hetero model
	@Test
	public void ringSplitHeteroAtom_rule4294_u78282()
	{
		String smirks = "[O:8]=[c:2]1[c:3][c:4][c:5][c:6][o:1]1>>[#8:1]([H])\\[#6:6]=[#6:5]/[#6:4]=[#6:3]\\[#6:2](-[#8-])=[O:8]";
		String smi = "C1=CC=C2C(=C1)C=CC(=O)O2";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	// fixed as of 2016-06-22 (3.0.3-SNAPSHOT/ambit2-smarts-3.0.3-20160621.122749-19.jar)
	@Test
	public void stereoChemLost_rule3138_u138720()
	{
		String smirks = "[#8:1]([H])-[#6:2](-[#6:9](-[#8-:10])=[O:11])=[#6:3](-[#1,#6,#17:12])-[#6:4]=[#6:5]-[#6](-[#8-])=O>>[#8-:10]-[#6:9](=[O:11])-[#6:2](=[O:1])-[#6:3](-[#1,#6,#17:12])-[#6:4]=[#6:5]";
		String smi = "C(=C(/C(=O)[O-])\\Cl)/C=C(\\C(=O)[O-])/O";
		String expectedSmiles;
		try
		{
			IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
					.parseSmiles("[O-]C(=O)C(=O)C\\C=C\\Cl");
			expectedSmiles = SmilesGenerator.absolute().create(mol);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		List<String> s = applySmirks(smirks, smi);
		Assert.assertEquals(expectedSmiles, s.get(0));
	}

	//apparently (!) fixed by nick
	@Test
	public void stereoChemIsLost2_rule1196_u45008()
	{
		String smirks = "[H][#6:1](-[#6:5])=[O:4]>>[#6:5]-[#6:1](-[#8-])=[O:4]";
		String smi = "C(=O)[C@@H]1C(=O)C(C(=O)O1)(F)F";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should still contain stereocheminfo: " + s.get(0),
				s.get(0).contains("@"));
		//		String smi2 = "C(=O)[C@@]([H])1C(=O)C(C(=O)O1)(F)F";
		//		List<String> s2 = applySmirks(smirks, smi2);
		//		Assert.assertTrue("Results should still contain stereocheminfo: " + s.get(0) + " " + s2.get(0), s.get(0)
		//				.contains("@") && s2.get(0).contains("@"));
	}

	// fixed by nick
	@Test
	public void invalidStereoChemProduct_rule4120_u134670()
	{
		String smirks = "[#8-:15]-[#6:1](=[O:16])\\[#6:2]([H])=[#6:3]([H])/[#6:4](=[#6:5]([H])\\[#6:6](-[#8-:8])=[O:7])/S([#8-])(=O)=O>>[#8-:15]-[#6:1](=[O:16])-[#6:2]-[#6:3]-[#6:4](-[#8-])=O.[#6:5]-[#6:6](-[#8-:8])=[O:7]";
		String smi = "C(=C/C(=O)[O-])/C(=C\\C(=O)[O-])/S(=O)(=O)[O-]";
		Assert.assertNotNull("Reaction should not fail", applySmirks(smirks, smi));
	}

	// as discussed with email in nick:
	// actually: implict Hs should not be added to products
	// instead: all smirks should state more explicitly what to do with Hs 
	@Test
	public void retainUnfilledValence_rule1196_u133400()
	{
		// the unfilled valence of nitrogen (-> [N]) should remaine unchanged
		String smirks = "[H][#6:1](-[#6:5])=[O:4]>>[#6:5]-[#6:1](-[#8-])=[O:4]";
		String smi = "C[N]C(=O)C(=O)C=O";
		List<String> s = applySmirks(smirks, smi);
		//Assert.assertTrue("Results should still contain [N]: " + s.get(0), s.get(0).contains("[N]"));

		// however, at the same time, Hs should be added to newly created atoms  
		smirks = "[H:5][C:1]([#6:6])([#1,#9,#17,#35,#53:4])[#9,#17,#35,#53]>>[H:5][C:1]([#6:6])([#8])[#1,#9,#17,#35,#53:4]";
		smi = "C(CN(CCCl)CC(C(=O)O)N)Cl";
		s = applySmirks(smirks, smi);
		Assert.assertFalse("Results should NOT contain [O]: " + s.get(0), s.get(0).contains("[O]"));
	}

	// fixed
	@Test
	public void stereoChemIsLost_rule4212_u136348()
	{
		String smirks = "[H:6][C:1]([#6:4])([#16;H1v2])[#1,#6:5]>>[H:6][C:1]([H])([#6:4])[#1,#6:5]";
		String smi = "CN\\C(NCCS)=C\\[N+]([O-])=O";
		List<String> s = applySmirks(smirks, smi);
		Assert.assertTrue("Results should still contain stereocheminfo: " + s.get(0),
				s.get(0).contains("\\"));
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
		Assert.assertFalse("Results should not contain C(=C=C: " + s.get(0),
				s.get(0).contains("C(=C=C"));
	}

	// fixed by rewriting smi(les): replaced C\\2\\C with C2\\C
	@Test
	public void productToSmilesError_rule3707_u143203()
	{
		String smirks = "[H:10][#8:9]-[c:4]1[c;R1:5][c;R1:6][c:1]([H])[c;R1:7][c;R1:8]1>>[H:10][#8:9]-[c:4]1[c;R1:5][c;R1:6][c:1](-[#8])[c;R1:7][c;R1:8]1";
		String smi = "C1=CC=C(C=C1)/C(=C2\\C=C(\\C(=C(\\C3=CC(=C(C=C3)O)O)/C(=O)[O-])\\C=C2O)O)/C(=O)[O-]";
		Assert.assertNotNull("Reaction should not fail", applySmirks(smirks, smi));
	}
	
	// fixed by replacing three single bonds with default (single or aromatic) bonds
	@Test
	public void missingProducts_rule2793_u26103()
	{
	// original smirks: "[#8:8]([H])-[c:2]1[c:1](-[#8:7]([H]))[c;R]([c;R:5](-[!#8,#1:11])[c;R:4](-[!#8,#1:10])[c;R:3]1-[!#8,#1:9])S([#8])(=O)=O>>[!#8,#1:11]\\\\\\[#6:5]=[#6:4](///[!#8,#1:10])-[#6:3](-[!#8,#1:9])-[#6:2](=[O:8])-[#6:1](-[#8-])=[O:7]";
		String smirks = "[#8:8]([H])-[c:2]1[c:1](-[#8:7]([H]))[c;R]([c;R:5]([!#8,#1:11])[c;R:4]([!#8,#1:10])[c;R:3]1[!#8,#1:9])S([#8])(=O)=O>>[!#8,#1:11]\\[#6:5]=[#6:4](/[!#8,#1:10])-[#6:3](-[!#8,#1:9])-[#6:2](=[O:8])-[#6:1](-[#8-])=[O:7]";
		String smi = "C1=C(C=CC(=C1)C2=NC3=C(C(=C(C=C3N2)S(=O)(=O)[O-])O)O)C4=NC5=C(C=C(C=C5S(=O)(=O)[O-])S(=O)(=O)[O-])N4";
		// working smiles example: C1=CC(=C(C(=C1N)S(=O)(=O)[O-])O)O
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	// fixed with ambit2-smarts-3.0.3-20160628.065959-20.jar
	@Test
	public void ringSplit1_alt_rule4224_arom_ar13()
	{
		String smirks = "[#8;H1:2]-[#6:10]1:[#6:5](-[#8;H1:1]):[#6:6](-[*,#1:11]):[#6:7]:[#6:8]:[#6,#7:9]:1>>[#8&-:2]-[#6:10](=O)-[#6,#7:9]=[#6:8]-[#6:7]=[#6:6](-[*,#1:11])-[#6:5](-[#8&-:1])=O";
		String smi = "Oc1cccc(Cl)c1O";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
	}

	// ambit2-smarts-3.0.3-20160628.065959-20.jar and using the aromatic smirks form while keeping the smiles kekulized
	@Test
	public void ringSplit1_alt_rule4224_kekulized_ar13()
	{
		String smirks = "[#8;H1:2]-[#6:10]1:[#6:5](-[#8;H1:1]):[#6:6](-[*,#1:11]):[#6:7]:[#6:8]:[#6,#7:9]:1>>[#8&-:2]-[#6:10](=O)-[#6,#7:9]=[#6:8]-[#6:7]=[#6:6](-[*,#1:11])-[#6:5](-[#8&-:1])=O";
		String smi = "OC1=CC=CC(Cl)=C1O";
		Assert.assertFalse("Results should not be empty", applySmirks(smirks, smi).isEmpty());
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
		String smirks = "[#6:4]@-&!:[#6;!$(C1(=O)C=CC(=O)C=C1)!$(C(=O)CC=O):1](@-&!:[#6:2])=[O:5]>>[#6:2]@-[#8]@-[#6:1](@-[#6:4])=[O:5]";
		// simple true target
		String smi = "O=C1CCOC=C1";
		List<String> s = applySmirks(smirks, smi);
		Map<String, String> diff = diffProductToReference(s, new String[]{"C1=COCCOC1=O","C1COC=COC1=O"});
		Assert.assertTrue(""+diff, null == diff);
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
		return TestAmbit.applySmirks(smrk, smi);
	}
	public static Map<String,String> diffProductToReference(List<String> products, String[] reference)
	{
		return TestAmbit.diffProductToReference(products, reference);
	}
}
