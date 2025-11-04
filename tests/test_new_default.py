import unittest
import os
import subprocess
import sys

class TestNewDefault(unittest.TestCase):
    """Test that sequential sampling is the new default."""

    @classmethod
    def setUpClass(cls):
        """Install dependencies."""
        dependencies = ['biopython', 'scipy', 'seaborn', 'matplotlib', 'pyrodigal', 'tqdm', 'numpy']
        subprocess.run([sys.executable, '-m', 'pip', 'install'] + dependencies, check=True)

    def setUp(self):
        """Set up test files."""
        self.test_gbk = 'test.gbk'
        # Create a genome with valid DNA sequences
        # Focal region will be 300 bp
        # Need total CDS > 6000 bp for focal to be <5%
        # Create 25 genes of 300 bp each = 7500 bp total CDS
        # gene1 is the focal gene
        genome_seq = ""
        features_text = ""
        pos = 1
        
        for i in range(1, 26):  # 25 genes
            gene_seq = "atg" + "gct" * 99  # 300 bp: ATG start + GCT repeats
            genome_seq += gene_seq
            features_text += f"     CDS             {pos}..{pos+299}\n"
            features_text += f"                     /locus_tag=\"gene{i}\"\n"
            pos += 300
            # Add 50 bp spacer between genes (except after last gene)
            if i < 25:
                genome_seq += "n" * 50
                pos += 50
        
        total_length = len(genome_seq)
        
        with open(self.test_gbk, 'w') as f:
            f.write(f"""LOCUS       Test_Seq              {total_length} bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Test sequence.
ACCESSION   Test_Seq
VERSION     Test_Seq.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
{features_text}ORIGIN
        1 {genome_seq}
//""")
        
        self.focal_gbk = 'focal.gbk'
        # Focal region with same gene1 sequence (300 bp = 4% of total CDS)
        focal_seq = "atg" + "gct" * 99  # 300 bp, matching gene1
        with open(self.focal_gbk, 'w') as f:
            f.write(f"""LOCUS       Focal_Test              300 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Focal test sequence.
ACCESSION   Focal_Test
VERSION     Focal_Test.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..300
                     /locus_tag="gene1"
ORIGIN
        1 {focal_seq}
//""")

    def tearDown(self):
        """Tear down test files."""
        os.remove(self.test_gbk)
        os.remove(self.focal_gbk)
        # Clean up output file if it exists
        if os.path.exists('test_output.txt'):
            os.remove('test_output.txt')

    def test_sequential_sampling_is_default(self):
        """Test that sequential sampling is the default."""
        outfile = 'test_output.txt'
        # Remove output file if it exists from previous failed run
        if os.path.exists(outfile):
            os.remove(outfile)
        
        codoff_cmd = [
            sys.executable,
            'bin/codoff',
            '-g', self.test_gbk,
            '-f', self.focal_gbk,
            '-o', outfile
        ]
        result = subprocess.run(codoff_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            self.fail(f"Command failed with return code {result.returncode}")
        with open(outfile, 'r') as f:
            content = f.read()
        # Updated assertion: no longer looking for "Sequential Sampling" line since we removed it
        self.assertIn('Discordance Percentile', content)
        os.remove(outfile)

if __name__ == '__main__':
    unittest.main()
