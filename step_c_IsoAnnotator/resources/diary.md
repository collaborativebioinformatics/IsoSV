### 2025/08/28

- Wrote a script called exon_check.py
- What does it does?
- It checks for the different SVs where could be its functional consequences
- Currently three funcions are implemented
	- if the sv lies within the predifiened exonic regions
	- if the sv lies completely within the exonic or partial in the exonic region
	- if it lies in the intergenic region 
	- if it is around the predefined splice site
	- Tries to get the final results as in the vcf format
- tried to get best mock data and generated vcf file based on it

Follow ups:
- Could be really interesting to work with the cleaned gtf file from the lui..
- The annotation can be really improved to find define tss and other functional regions
- Codes could be cleaned up and typed for better organization and reading

Louis progress: 
- Fixed a bug in tx tree construction, tx_tree can now pinpoint exact overlapping gene regions
- Added code based on the reconstructed tx tree for SV annotation (INS associated gene fusions will need extensive work into BAM?)
	- Current code in `annotate.py`, still in tsv format
- Not yet converted to VCF format
- Not yet tested on most up-to-date step b files
- Refactored the old 01 02 scripts 