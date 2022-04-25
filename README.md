## MetBP : *A standalone command line tool for detection and analysis of metal-ion and Basepair interactions*


## Synopsis  


> Metal Basepair interaction program is a software tool that detects and analyzes metal interactions with RNA (It can handle DNA anso) base pairs from the crystal structure files stored in mmCIF or PDB format. The program gives results in suitable format in plain text as well as in machine readible JSON and CSV formats also. The program is a stand alone command line based tool developped for Linux. The program is written in C and FORTRAN.

## Installation 
>Download at least all the files of sys and bin directory.
	Transfer the 'bin/metallic' executable file to some path folder
	and transfer all the files of *'sys'* folder to some suitable location.
	We prefer the following.

>		sudo cp bin/metbp.linux   /usr/local/bin/
	
>		sudo cp sys/*   /usr/local/bin/


## Compilation
> run make. The executable file will be stored in the bin directory.
## Dependencies
> The compilation process needs a C and a FORTRAN compiler. We have assumed *gcc* for C and *gfortran* for FORTRAN. We have made the Makefile accordingly. If the user wants a different compiler, change the makefile accordingly. 

## About The Makefile
> 
We have given a ready made make file for easy compilation.  The Makefile assumes *gcc* as C compiler and *gfortran* as FORTRAN compiler. If, however, the user do not have access to these compilers and to something else instead, they need a small modification in the Makefile. For example, if you have *clang* as C compiler and *f77* as fortran compiler, then chage the following in the Makefile. 
	change `CC := gcc -std=c99` to `CC := clang` and change `FF := gfortran` to `FF := f77`

##  Basepairs
>Nucleic acids forms base-base interactions for their stebelization. These base-base interactions are called basepairs. In DNA we mainly observe that a Guanine forms a basepair with a Cytosine and an Adenine forms a baspair with Thymie. In RNA, however, varieties of basepairs are observed. Typically 158 types of base pairs are theoratically possible, though only 126 of them are found in nature. Out of them five base pairs are called canonical and the rest are called non-canonical. The canonical basepairs are GC cWW, AU cWW and GU   cWW. 

##  A note about basepair writing convension
>The most popular notation for basepair is due to Leontis-Westhof. Their notation looks like this.
>
		     GC cHW (Leontis-Westhof nomenclature)
This means a Guaning forms a basepair with Cytosine in Cis orientation and Guaninen’s Hoogsteen edge and Cytosin’s Watson-Crick edge forms the basepair. In our case, we write the the same as follows,
>
	         G:C H:WC (BPFind nomenclature that is followed by MetBP)
> The reason is in MetBP program, we consider C-H…O/N mediated basepairs as separate basepairs and we denote them in small letters. So, in MetBP G:C H:WC and G:C h:wC carries different meaning. G:C H:WC is a normal basepair whereas G:C h:wC is a basepair that contains C-H…O/N interations. Leontis-Westhof nomenclature does not differentiate them. Futrher, MetBP program considers protonated base pairs also and reports them differently. Following is the full base pair nomenclature that MetBP follows.
>	
			W - Watson-Crick edge (Capital W).
			H - Hoogsteen edge (Capital H).
			S - Sugar edge (Capital S).
			w - Watson-Crick edge with one C-H...O/N type of hydrogen bond (Small w).
			h - Hoogsteen edge with one C-H...O/N type of hydrogen bond (Small h).
			s - Sugar edge with one C-H...O/N type of hydrogen bond (Small s).
			+ - Protonated Watson-Crick edge.
			z - Protonated Sugar edge.
			g - Protonated Hoogsteen edge (very rare in nature).

## Basepair types
> BPNet used three types of basepair types. These are as follows
>
	BP - Normal base pair.
	TP - Tartiary pair.
	BF - Bifurcated pair indicating a single edge of the central base is paired to two bases simultaneously 
	
([Following our paper BPFIND,2006](https://doi.org/10.1080/07391102.2006.10507108))
	

## Basepair rule
In MetBP, we follow the BPFIND algo


## Command-line options

>###  Software Related command-line options
>>**--help:** Shows a small help on the terminal.<br> 
**--genhelp:** Generates this help.md file in the current directory. <br>
**--version:** Prints the version of the program.<br>
**--pub:** Shows how to refer to the software if it is used by anyone.<br>
**--contact:** Prints the contact information for query or bug report.<br>


> ### Domain related command line options.
>> 
**–hbdist=[value]:** MetBP uses 3.8A as default distance for donor-acceptor atom distance for hydrogen bonding interactions. However, this can be changed with this option.  <br><br>
**EXAMPLE** <br>
    `metbp.linux 1n32.cif -hbdist=4.1`<br><br><br>
**-sugmed=[true/false]:** By default MetBP considers pentose sugar O2’ atom for sugar edge based pairing. This is called a sugar mediated base pair. If the user does not want that, then they can set it to ‘false’. <br><br>
**EXAMPLE**<br>
`metbp.linux 1n32.cif -sugmed=false` <br><br><br>
**-chmed=[true/false]:** By default MetBP considers C-H...O/N mediated hydrogen bonds for base pair formations. If the user does not want it, the C-H…O/N type hydrogen bonds can be excluded. <br><br>
**EXAMPLE** <br>
`metbp.linux 1n32.cif -chmed=false`<br><br><br>
**-hetatm=[true/false]:** By default, the MetBP software considers all modified nucleic acids. The modified nucleic acid residues can be excluded from computations with the following command line option. <br><br>
**EXAMPLE**<br>
`metbp.linux 1n32.cif -hetatm=false`<br><br><br>
**-chain=name:** The MetBP program can handle a specific chain also. If the user supplies the chain name, the program considers base pairs of only that specific chain and thus metal and base pair binding of that specific chain.<br><br>
**EXAMPLE** <br>
`metbp.linux 1n32.cif -chain=A`<br><br>
**CAUTION !!!**  *If a specific chain is selected then all analyses are restricted to that specific chain only. So, if a base of the given chain has a pair with the bases of the other chain, the program will not report them.* <br><br><br>
**-nmrmdl=[number]:** If an NMR structure is supplied, by default the MetBP software reads the first model and works on that. But it can take a specific model of an NMR structure and can compute the metal and base pair interaction on that model also. <br><br><br>
**EXAMPLE** <br>
`metbp.linux 2ll9.cif -nmrmdl=4`<br><br><br>
**-mode=[bp/nuc/all/dev]:** The MetBP program works in three different modes. The ‘bp’ mode, the ‘nuc’ mode and finally the ‘all’ mode.  In ‘bp’ mode, the details of only those metals are reported which are directly coordinating with an atom of nucleic acid (i.e. base part) that forms a base pair with some other bases. (Note: It may be noted that we consider O2’ of pentose sugar as an important atom for base pairs. So, any bond with O2’ will be treated as an atom of base.) In ‘nuc’ mode, it reports the details of all metals which bind with any atom of a nucleic acid residue. i.e. in this mode, if a metal binds with the oxygen of a phosphate group or any atom of the pentose sugar, the details of that metal is considered. In 'dev' mode, the program runs internally in 'bp' mode but it generates two extra json files. One contains all base pair information and the other contains all metals that has a contact with nucleic acids (The nucleic acis may not form a base pair). <br><br><br>
**-paramfile=[PATH]:** Default the program takes the metal parameters from the ‘metal.params’ file stored in the NUCLEIC_ACID_DIR path location. If the user wants to alter that metal file, she/he may do that. However, if the user does not want to modify the default metal.params file, she may take a copy of the same in any other location and can supply the same by this switch. <br><br>
**EXAMPLE** <br>
`metbp.linux 1n32.cif -paramfile=./metal.params`<br><br>
**Note:**  *If the user supplies the parameter file in this way, then all the parameters will be taken from this file. * <br><br>
**Caution!!!** *This is not an update of the existing values. It is the resetting all the values. For example, if the user supplies the values for Magnesium, then the system will consider only Magnesium and no other metals.* <br><br><br>



| Metal Symbol | Atomic No | Oxygen Dist | Nitrogen Dist | H2O Dist | Carbon Dist | Sulfur Dist | Phosp Dist | Metal Dist|
|---------|---------|---------|---------|---------|---------|---------|---------|---------
|   Na   | 11      | 2.90  | 3.1    | 3        | ?        |?         | ? | 2.9|

























##  List of output files

- **1n32.sum:** This is the summary file. This file reports the number of nucleic acids and protein residues, metals, water molecules found in the structure. Also it gives a summary of every metal binding site, i.e. its coordination, number of nucleic acids, proteins etc to which it binds. <br><br>
- **1n32.met:** This is the main output file. This file gives the metal and base pair interaction details. This file also gives the secondary structure of the base with which the metal binds.<br><br>
- **1n32.det:** This file gives the details of every metal binding site. It gives the metal to the coordinated atom distances, different atom-metal-atom distances.<br><br>
- **1n32.pml:** MetBP generates a [pymol](https://pymol.org/) script. The script, when run through pymol, shows the metal and its different binding residues. To give a better visual effect, the different types of bases are shown in different colors. For example, white bases indicate that they form a base pair. Orange colored base indicates that the  metal binds to any atom of the pentose sugar. Green colored base indicates that the metal binds to the oxygen of the phosphate group.<br><br>
- **1n32_metbp.json:** This file gives the metal-base pair interaction details in JSON format.<br><br>
 - **1n32_basepair.json:** This file gives all the base pair details in JSON format. The get this file the user has to run the program in developer mode . i.e. the program has to be run as
 `metbp.linux 1n32.cif  -mode=dev`<br><br>
 - **1n32_metnuc.json:** This file gives the information of all metals that interact with a nucleic acid resisue whether or not that residue forms a pair. The get this file the user has to run the program in developer mode . i.e. the program has to be run as
 `metbp.linux 1n32.cif  -mode=dev`<br><br>

## run:
>
`metbp.linux [OPTIONS] [mmCIF/PDB]`<br>
**EXAMPLE:**<br>
`metbp.linux 1ehz.cif`<br>

## bug-report:

		Parthajit Roy, 
			`roy.parthajit@gmail.com`

		Dhananjay Bhattacharyya,
			`dhananjay.bhattacharyya@saha.ac.in` 
