# MGFdissect  
## Short description  
*"dissecting" and reordering MGF*
   
The tool in this main branch parses and groups MGF (mascot generic format) MS/MS spectra from multiple files into a single file. The reordering and grouping of particular entries (between BEGIN IONS and END IONS statements) is based on their parent ion m/z and retention time values. All spectra are kept in the merged file, except for 1- and 0-liners. A minimalistic user interface is provided through svDialogs package.  
  
## Usage  
  
Mass spectrometry instrument vendor software can export single or multiple files as .mgf format. However, the focus is on samples rather than individual features (compounds). Whenever dealing with a series of related samples with the goal to identify compounds present in one or multiple samples it is useful to group the features together (same compounds appear in multiple samples), even before the identification step. This is exactly what MGFdissect does. This kind of data pretreatment can speed up the evaluation/identification with MS/MS libraries immensely.   
  
## Workflow  
  
After exporting data dependent MS/MS measured chromatograms (e.g. AutoMS) to .MGF format the script (mgf-diss.R) is run while following the prompts. The generated new merged mgf is automatically saved in the same folder as the original series of .mgf files (as "merged.mgf"). Example files are in the mgfs folder.  
  
You can learn more about the tool in the Wiki 