import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))
sys.path.append(str(os.path.dirname(__file__)) + "/../../workflow/scripts")

import common

def test_get_ncbi_sequences():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/integration/get_ncbi_sequences/data")
        expected_path = PurePosixPath(".tests/integration/get_ncbi_sequences/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # write run_assembly.txt in temporary workdir
        with open(str(workdir) + "/run_assembly.txt", "w") as run_assembly:
            run_assembly.write("assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material	asm_not_live_date\n")
            run_assembly.write("GCF_001735525.1	PRJNA224116	SAMN05384437	MCBT00000000.1	representative genome	23	23	Shewanella colwelliana	strain=CSB03KR		latest	Scaffold	Major	Full	2016/09/19	ASM173552v1	Chonnam National University	GCA_001735525.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/735/525/GCF_001735525.1_ASM173552v1			na")

        # Run the test job
        import downloadGenes as dg
        os.system("mkdir " + str(workdir) + "/ncbi/")
        downloaded_genome_ids, failed_genome_ids = dg.download_genomes(str(workdir))

        if failed_genome_ids:
            raise ValueError("Failed to download: " + str(failed_genome_ids))
        for id in downloaded_genome_ids:
            if (not os.path.isfile(str(workdir) + "/ncbi/" + id + "_cds_from_genomic.fasta")) or (not os.path.isfile(str(workdir) + "/ncbi/" + id + "_rna_from_genomic.fasta")):
                raise ValueError("Couldn't find: " + id)