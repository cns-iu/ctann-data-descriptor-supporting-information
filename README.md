# Constructing, Validating, and Using Human Reference Atlas Cell Type Annotation Crosswalks

Nicole Vasilevsky1,2*, Aleix Puig-Barbe3*, Andreas Bueckle1, Ellen M. Quardokus1, Yongxin Kong1, Jie Zheng4, Yashvardhan Jain1, Bruce W. Herr II1, Katy Börner1*

1 Department of Intelligent Systems Engineering, Luddy School of Informatics, Computing, and Engineering, Indiana University, Bloomington, IN, USA
2 Rose City Data Science, Portland, OR, USA
3 European Bioinformatics Institute (EMBL-EBI), Wellcome Genome Campus, Hinxton, Cambridge CB10 1SD, UK
4 University of Michigan Medical School, Ann Arbor, MI, USA

## Abstract

The Human Reference Atlas (HRA) v2.5 includes 4,807 anatomical structures linked to Uberon or FMA, 2,215 cell types, and 1,347 biomarkers (genes, proteins, lipids); 73 organs with 1,356 anatomical structures are modeled in three dimensions (3D). Mapping new data into the HRA requires 3D spatial registration and/or cell type annotation (CTann). Many CTann tools exist but few use cell type names crosswalked to CL—making comparisons across tools, tissues, and organs difficult. The HRA supports five transcriptomics CTann tools (Azimuth, CellTypist, popV, FR-Match, Pan-Human Azimuth), four proteomics CTann tools (DeepCell Types, DeepCell Types-HuBMAP, Robust Image-Based Cell Annotator [RIBCA], and SpaTial cELl LeARning [STELLAR]), plus a CDE Spatial Omics collection covering previously published annotations across 12 tissue types with a total of 47,349,496 cells. Assigned cell types differ in number, resolution, and coverage. This paper details how 10 crosswalks were constructed,  validated, and used to compile a ‘Supertree’ for cross-organ and technology CTann. We conclude with recommendations for CTann developers that will ease crosswalk generation in the future. All data and code is freely available at https://github.com/cns-iu/ctann-data-descriptor-supporting-information.  
