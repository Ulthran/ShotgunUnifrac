import os
import tempfile

TESTS_FP = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_FP = os.path.join(TESTS_FP, "test-data")
TEMP_FP = tempfile.mkdtemp()
OUTPUT_FP = os.path.join(TEMP_FP, "output/")
GENOMES_FP = os.path.join(OUTPUT_FP, "genomes/")
FILTERED_FP = os.path.join(OUTPUT_FP, "filtered-sequences/")
FILTERED_NUCL_FP = os.path.join(OUTPUT_FP, "filtered-nucl-sequences/")
MERGED_FP = os.path.join(OUTPUT_FP, "merged-sequences/")
