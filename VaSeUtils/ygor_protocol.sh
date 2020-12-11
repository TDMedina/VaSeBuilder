#string ygorLocation
#string vbGoldStandard
#string vbVcfCalls
#string ygorOutdir
#string project

echo "Compare VaSeBuilder variant calls to gold standard variant calls"
python "${ygorLocation}" -1 "${vbGoldStandard}" -2 "${$vbVcfCalls}" -o "${ygorOutdir}" -op "${project}"
echo "Comparison to gold standard complete"

#SHORT EXPLANATION ON THE DEFINED PARAMETERS
#ygorLocation = Path to the location of ygor.py tha will perform the comparison
#vbGoldStandard = Path to VCF file containing the original calls of all variants added to the VaSeBuilder sample
#vbVcfCalls = Path to the VCF file for calls on the VaSeBuilder sample produced by the NGS_DNA pipeline
#ygorOutDir = Path to the directory where the comparison output files should be written to
#project = The project name used by the pipeline. This name will be used in the name of the comparison output files
