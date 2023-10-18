# MGFdissect  
## Short description  
*"dissecting" and reordering MGF*
   
The tool in this branch parses and groups MGF (mascot generic format) MS/MS spectra from multiple files into a single file. The reordering and grouping of particular entries (between BEGIN IONS and END IONS statements) is based on their parent ion m/z and retention time values. All spectra are kept in the merged file, except for 1- and 0-liners. A minimalistic user interface is provided through svDialogs package.     
This script apart from combining the MFGs searches also a custom .csv database with a similar structure as [mzmine](https://github.com/mzmine/mzmine_documentation) custom database (example included in the repo). 