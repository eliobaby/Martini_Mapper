def get_m3_dict():
    # 1) Start by listing all (key, value) pairs in one big list.
    #    Even if a key repeats in this list, we’ll handle collisions when inserting into the dict.
    entries = [
        # -------------------------
        # Section 1: Benzene nodes
        # -------------------------
        
        # FROM THE BUILDING BLOCK TABLE: (unchangable)
        ("TN6a", [1, 0, "CN", 2, "", []]),
        ("TC5",  [1, 0, "CC", 2, "", []]),
        # FROM TABLE 24: (no change)
        # aromatic (no hydrogens)
        # ("TC5e",  [1, 0, "CC", 2, "", []]),
        
        # FROM 90 ORIGINAL MOLECULES:
        ("TN1a",  [1, 0, "NN", 2, "", []]), #TN3a

        # ----------------------------------
        # Section 2: Benzene border nodes
        # ----------------------------------
        
        # FROM THE BUILDING BLOCK TABLE: (unchangable)
        ("SN6d",  [2, 0, "CC(N)",   2, "", []]), #SN4d 
        ("TN6a",  [2, 0, "C(=O)",   2, "", []]), #TN2a 
        ("SN2a",  [2, 0, "C(OC)",   2, "", []]), #SN1a 
        ("SC6",   [2, 0, "CC(S)",   2, "", []]),
        ("SC4",   [2, 0, "CC(C)",   2, "", []]),
        ("TC4",   [2, 0, "C(C)",    2, "", []]), #TC6 
        ("SX3",   [2, 0, "CC(Cl)",  2, "", []]),
        ("SX2",   [2, 0, "CC(Br)",  2, "", []]), #SX1 
        ("X1",    [2, 0, "CC(I)",   2, "", []]),
        
        # FROM TABLE 24: (no change)
        # Phenol
        ("SN6",   [2, 0, "CC(O)",   2, "", []]), #SN2 after old tuning
        ("TN6+",  [2, 0, "C(O)",    2, "", []]), #TN2 
        # Phenol + Phenol = Diol
        ("SP4",   [2, 0, "C(O)(CO)", 2, "", []]),
        # FROM 90 ORIGINAL MOLECULES:
        # methyl pyrrole
        ("TN1",   [2, 0, "N(C)",  2, "", []]), #TN3 
        
        # INFERRED AND TUNED from old model
        ("SX4e",  [2, 0, "CC(F)",   2, "", []]), #SX3e
        ("SN6a",  [2, 0, "CC(=O)",  2, "", []]), #SP1a
        ("TX3",   [2, 0, "C(Cl)",   2, "", []]), #TX4
        ("TN6d",  [2, 0, "C(N)",    2, "", []]),
        ("TC6",   [2, 0, "C(S)",    2, "", []]),
        ("TX1",   [2, 0, "C(I)",    2, "", []]),
        ("TX2",   [2, 0, "C(Br)",   2, "", []]), #TX1
        ("TX4e",  [2, 0, "C(F)",    2, "", []]),
        
        
        # ----------------------------------
        # Section 3: Non‑benzene 6‑ring nodes
        # ----------------------------------
        
        # FROM THE BUILDING BLOCK TABLE: (unchangable)
        ("SN4a", [3, 0, "COC", 2, "SC3"]),
        ("SN3a", [3, 0, "COC", 2, "SN3a"]),
        ("SC3",  [3, 0, "CCC", 2, "", []]),
        
        # FROM TABLE 24: (no change)
        # primary amine
        ("SN6d", [3, 0, "CCN", 2, "", []]),
        ("TN6d", [3, 0, "CN", 2, "", []]),
        
        # -----------------------------------------        
        # Section 4: Non‑benzene 6‑ring border nodes
        # -----------------------------------------
        # INFERRED AND TUNED from old model
        ("SN5a", [4, 0, "OCO", 2, "", []]),
        ("SN1",  [4, 0, "NCN", 2, "", []]),
        ("SX4e", [4, 0, "FCF", 2, "", []]),
        ("SN3a", [4, 0, "COC", 2, "", []]),
        # FROM 90 ORIGINAL MOLECULES:
        # methyl pyrrole
        ("TN1",  [4, 0, "N(C)",  2, "", []]), #TN3 
        # -----------------------
        # Section 5: 5‑ring nodes
        # -----------------------
        
        # FROM THE BUILDING BLOCK TABLE: (unchangable)
        ("TN6d", [5, 0, "NH",  2, "", []]),
        ("SN5a", [5, 0, "OCO", 2, "", []]),
        ("SC3",  [5, 0, "CCC", 2, "", []]),
        ("TN4a", [5, 0, "OC",  2, "", []]), #TP6a
        ("TN2a", [5, 0, "O",   2, "", []]), #TN3a 
        ("TC3",  [5, 0, "CC",  2, "", []]),
        # Thiophene
        ("SC6",  [5, 0, "CSC", 2, "", []]),
        ("TC6",  [5, 0, "S",   2, "", []]), #TN2 
        
        # INFERRED AND TUNED from old model
        ("TC6",  [5, 0, "SS",  2, "", []]),

        # ------------------------------
        # Section 6: 5‑ring border nodes
        # ------------------------------
        
        # FROM THE BUILDING BLOCK TABLE:(unchangable)
        ("SC3", [6, 0, "CC(C)", 2, "", []]),
        # methyl pyrrole
        ("TN1", [6, 0, "N(C)",  2, "", []]), #TN3 
        
        # FROM 90 ORIGINAL MOLECULES:
        ("SN6",  [6, 0, "CC(O)", 2, "", []]),
        ("TN6a", [6, 0, "CC(=O)", 2, "", []]),
        
        # INFERRED AND TUNED from old model
        ("SN6d", [6, 0, "CC(N)", 2, "", []]), #SP6
        
        # --------------------------------------------------
        # Section 7: Non‑ring nodes (polar/apolar fragments)
        # --------------------------------------------------
        
        # FROM THE BUILDING BLOCK TABLE: (unchangable)
        # tri-halogen-o-methane
        ("SX4e", [7, 1, "C(F)(F)(F)",      1, "", []]),
        ("X2",   [7, 1, "C(Cl)(Cl)(Cl)",   0, "", []]),
        # nitrile
        ("TN4a", [7, 0, "N#C",          1, "", []]),
        # Nitrite
        ("SN3a", [7, 0, "ON=O",         1, "", []]),
        # alcohol or ether
        # I put dummy data for val[6] bc conflicts
        ("TP1d", [7, 0, "CO",           1, "SX4e", ["O", 2]]), 
        ("TN2a", [7, 0, "CO",           1, "TC5", ["O", 2]]), #TN4a 
        # tertiary amine
        ("SN1",  [7, 0, "CNC",          1, "", []]), #SN1 
        # linear alkane
        ("TC3",  [7, 0, "CC",           1, "", []]),
        ("SC2",  [7, 0, "CCC",          1, "", []]), #SC3 
        # aldehyde
        ("N6a",  [7, 0, "CCC=O",        1, "C1"]), 
        ("TN4a", [7, 0, "C=O",          1, "", []]),
        ("SN4a", [7, 0, "CC=O",         1, "", []]), 
        #("SN3r", [7, 0, "COC",          0, "", []]), describe in ether
        #("P1",   [7, 1, "C(O)(C)(C)",     0, "", []]), describein alcohol
        #("N4a",  [7, 0, "COC=O",        1, "", []]), #N5a describe in ester
        #("C1",   [7, 0, "CCCC",         0, "", []]), #C3 describe in linear alkane
        #("C4",   [7, 0, "CC=CC",        0, "", []]), describe in alkene
        #("TC5",  [7, 0, "C=C",          1, "", []]), describe in dienes
        # describe in carboxylic acid
        #("P2",   [7, 1, "C(O)(=O)(CC)",   0, "", []]), #P3
        #("SP2",  [7, 1, "C(O)(=O)(C)",    0, "", []]), #SN5
        # describe in ester
        #("N4a",  [7, 1, "C(C)(=O)(OC)",    0, "", []]), #N5a
        #("N4a",  [7, 0, "COC=O",        1, "", []]), #N5a
        
        # FROM 90 ORIGINAL MOLECULES:
        ("SC1",  [7, 1, "C(C)(C)(C)",     0, "", []]), # would have been SC2 by table 24 standard
            
        # FROM GRUNEWALD'S DATASET:
        # Sulfonamide
        ("P4",   [7, 1, "S(=O)(=O)(N)",   0, "", []]),
        # N-oxy sulfonamide
        ("P4",   [7, 1, "S(=O)(=O)(NO)",   0, "", []]),
        # N-carbon sulfonamide
        ("P4",   [7, 1, "S(=O)(=O)(NC)",   0, "", []]),
        # Cyanamide
        ("P4d",  [7, 0, "NC#N",           1, "", []]),
        # bromoalkene
        ("SX1",  [7, 0, "BrC=C",          1, "", []]),
        ("SX1",  [7, 0, "C=CBr",          1, "", []]),
        # aldehyde
        ("N4a",  [7, 1, "C(C)(C)(C=O)",   1, "", []]),
        # Sulfonyl fluoride
        ("X4",   [7, 1, "S(=O)(=O)(F)",   0, "", []]),
        # Sulfone
        ("SP5",  [7, 0, "O=S=O",         0, "", []]), #SP6?
        ("P4",   [7, 1, "S(=O)(=O)(C)",   0, "", []]), #SP3?
        # carbon disulfide
        ("C6",   [7, 0, "S=C=S",         0, "", []]), #SC6?
        # carbon dioxide
        ("SC6",  [7, 0, "O=C=O",         0, "", []]),
        # chlorobromomethane
        ("SX2",  [7, 0, "ClCBr",         0, "", []]),
        ("SX2",  [7, 0, "BrCCl",         0, "", []]),
        # chlorofluoromethane
        ("SX4",  [7, 0, "FCCl",          0, "", []]),
        ("SX4",  [7, 0, "ClCF",          0, "", []]),
        # diiodomethane
        ("X2",   [7, 0, "ICI",           0, "", []]),
        # Methylhydrazine
        ("P5d",  [7, 0, "CNN",           0, "", []]),
        # trichloroethylene
        ("X3h",  [7, 1, "C(=CCl)(Cl)(Cl)",0, "", []]),
        # vinylchloride
        ("SX3",  [7, 0, "C=CCl",         0, "", []]), #SX3
        ("SX3",  [7, 0, "ClC=C",         0, "", []]), #SX3
        
        # FROM TABLE 24: (no change)
        # "#" are old tuned version of the bead
        # linear alkane
        ("C1",   [7, 0, "CCCC",         0, "", []]),
        # branched alkane
        ("C2",   [7, 2, "C(C)(C)(C)(C)" , 0, "", []]), #C1
        #("SC2",  [7, 1, "C(C)(C)(C)",     0, "", []]), #SC3
        ("C2",   [7, 1, "C(CC)(C)(C)",    0, "", []]),
        # alkene
        ("SC4",  [7, 1, "C(=C)(C)(C)",    0, "", []]),
        ("C4",   [7, 1, "C(=CC)(C)(C)",   0, "", []]),
        ("C4",   [7, 1, "C(=C)(CC)(C)",   0, "", []]),
        ("C4",   [7, 1, "C(=C)(C=C)(C)",  0, "", []]),
        ("C4",   [7, 0, "C=CCC",        0, "", []]), #C6
        ("SC4",  [7, 0, "C=CC",         0, "", []]),
        ("C4",   [7, 0, "CC=CC",        0, "", []]),
        # dienes
        ("SC5",  [7, 0, "C=C=C",        0, "", []]), #SC4
        ("SC5",  [7, 0, "C=CC=C",       0, "", []]), #new
        ("TC5",  [7, 0, "C=C",          1, "", []]),
        # alkynes
        ("C6r",  [7, 0, "CC#CC",        0, "", []]),
        ("C6r",  [7, 0, "C#CCC",        0, "", []]),
        ("C6r",  [7, 0, "C#CC=C",       0, "", []]),
        ("C6r",  [7, 0, "C#CC#C",       0, "", []]),
        ("SC6r", [7, 0, "C#CC",         0, "", []]),
        ("TC6r", [7, 0, "C#C",          0, "", []]),
        # thiol/sulfide
        ("SC6",  [7, 0, "CCS",          2, "", []]),
        ("SC6",  [7, 0, "CSC",          0, "", []]),
        ("TC6",  [7, 0, "CS",           2, "", []]),
        ("SC6",  [7, 1, "C(C)(C)(S)",   0, "", []]),
        ("C6",   [7, 1, "C(C)(C)(SC)",  0, "", []]),
        ("C6",   [7, 1, "C(C)(C)(CS)",  0, "", []]),
        ("C6",   [7, 1, "C(C)(CC)(S)",  0, "", []]),
        # primary imine
        ("SN2",  [7, 0, "CC=N",         0, "", ["N", 1]]),
        ("TN2",  [7, 0, "C=N",          0, "", ["N", 1]]),
        ("SN2",  [7, 0, "C=C=N",        0, "", ["N", 1]]),
        ("SN2",  [7, 1, "C(C)(C)(=N)",  0, "", ["N", 1]]),
        # secondary imine
        ("SN1a", [7, 0, "CC=N",         0, "", ["N", 0]]),
        ("TN1a", [7, 0, "C=N",          0, "", ["N", 0]]),
        ("SN1a", [7, 0, "C=C=N",        0, "", ["N", 0]]),
        ("SN1a", [7, 1, "C(C)(C)(=N)",  0, "", ["N", 0]]),
        ("SN1a", [7, 0, "C=NC",         0, "", []]),
        ("N1a",  [7, 0, "CC=NC",        0, "", []]),
        ("N1a",  [7, 0, "C=NCC",        0, "", []]),
        # acyl chloride
        ("SN2ah",[7, 0, "O=CCl",        0, "", []]),
        ("SN2ah",[7, 0, "ClC=O",        0, "", []]),
        ("SN2ah",[7, 1, "C(C)(=O)(Cl)", 0, "", []]),
        ("N2ah", [7, 1, "C(CC)(=O)(Cl)",0, "", []]),
        # nitrile
        ("SN2a", [7, 0, "N#CC",         1, "", []]),
        # ether:
        ("SN3r", [7, 0, "COC",          0, "", []]),
        ("N3r",  [7, 0, "COCC",         0, "", []]),
        ("N3r",  [7, 0, "CCCO",         1, "", ["O", 0]]), #N2
        ("SN3r", [7, 0, "CCO",          1, "", ["O", 0]]), #SN2
        ("TN3r", [7, 0, "CO",           1, "", ["O", 0]]),
        ("N3r",  [7, 1, "C(O)(C)(CC)",  0, "", ["O", 0]]),
        ("N3r",  [7, 1, "C(C)(C)(CO)",  0, "", ["O", 0]]),
        ("N3r",  [7, 1, "C(O)(C)(C)",   0, "", ["O", 0]]),
        ("N3r",  [7, 1, "C(C)(C)(OC)",  0, "", []]),
        # alcohol
        ("P1",   [7, 1, "C(O)(C)(C)",   0, "", ["O", 1]]),
        ("P1",   [7, 1, "C(O)(C)(CC)",  0, "", ["O", 1]]),
        ("P1",   [7, 1, "C(C)(C)(CO)",  0, "", ["O", 1]]),
        ("P1",   [7, 0, "CCCO",         1, "", ["O", 1]]), 
        ("SP1",  [7, 0, "CCO",          1, "", ["O", 1]]), 
        ("TP1",  [7, 0, "CO",           1, "", ["O", 1]]), 
        # Ketone/aldehyde
        ("SN5a", [7, 1, "C(=O)(C)(C)",    0, "", []]),
        ("N5a",  [7, 1, "C(=O)(C)(CC)",   0, "", []]),
        ("N5a",  [7, 1, "C(=C)(C)(C=O)",   1, "", []]),
        ("N5a",  [7, 1, "C(C=C)(C)(=O)",   1, "", []]),
        # primary amin
        ("TN6d", [7, 0, "CN",           0, "", ["N", 2]]), #TP1
        ("SN6d", [7, 0, "CCN",          0, "", ["N", 2]]), #SN5
        ("SN6d", [7, 0, "C=CN",         0, "", ["N", 2]]), #SP1
        ("N6d",  [7, 0, "CCCN",         0, "", ["N", 2]]), #P3
        ("SN6d", [7, 1, "C(N)(C)(C)",   0, "", ["N", 2]]),
        ("N6d",  [7, 1, "C(N)(CC)(C)",  0, "", ["N", 2]]),
        ("N6d",  [7, 1, "C(N)(C=C)(C)", 0, "", ["N", 2]]),
        # secondary amine
        ("TN5",  [7, 0, "CN",           0, "", ["N", 1]]), 
        ("SN5",  [7, 0, "CCN",          0, "", ["N", 1]]), 
        ("SN5",  [7, 0, "C=CN",         0, "", ["N", 1]]), 
        ("N5",   [7, 0, "CCCN",         0, "", ["N", 1]]), 
        ("N5",   [7, 0, "CNCC",         0, "", ["N", 1]]), 
        ("N5",   [7, 0, "CNC=C",        0, "", ["N", 1]]), 
        ("SN5",  [7, 1, "C(N)(C)(C)",   0, "", ["N", 1]]),
        ("N5",   [7, 1, "C(N)(CC)(C)",  0, "", ["N", 1]]),
        ("N5",   [7, 1, "C(N)(C=C)(C)", 0, "", ["N", 1]]),
        ("N5",   [7, 1, "C(NC)(C)(C)",  0, "", ["N", 1]]),
        # tertiary amine
        ("TN3a",  [7, 0, "CN",           0, "", ["N", 0]]), 
        ("SN3a",  [7, 0, "CCN",          0, "", ["N", 0]]), 
        ("SN3a",  [7, 0, "C=CN",         0, "", ["N", 0]]), 
        ("N3a",   [7, 0, "CCCN",         0, "", ["N", 0]]), 
        ("N3a",   [7, 0, "CNCC",         0, "", ["N", 0]]), 
        ("N3a",   [7, 0, "CNC=C",        0, "", ["N", 0]]), 
        ("SN3a",  [7, 1, "C(N)(C)(C)",   0, "", ["N", 0]]),
        ("N3a",   [7, 1, "C(N)(CC)(C)",  0, "", ["N", 0]]),
        ("N3a",   [7, 1, "C(N)(C=C)(C)", 0, "", ["N", 0]]),
        ("N3a",   [7, 1, "C(NC)(C)(C)",  0, "", ["N", 0]]),
        ("SN3a", [7, 1, "N(C)(C)(C)",    0, "",  []]), #SN2a
        ("N3a",  [7, 1, "N(C)(C)(C=C)",  0, "", []]),
        ("N3a",  [7, 1, "N(C)(C)(CC)",   0, "", []]),
        # enol
        ("SP2",  [7, 0, "C=CO",          0, "", []]),
        ("P2",   [7, 0, "CC=CO",         0, "", []]),
        ("SP2",  [7, 1, "C(O)(=C)(C)",   0, "", []]),
        ("P2",   [7, 1, "C(O)(=CC)(C)",  0, "", []]),
        # carboxylic acid
        ("SP2",  [7, 0, "OC=O",         1, "", ["O", 0, 1]]),
        ("P2",   [7, 1, "C(O)(=O)(CC)", 0, "", ["O", 0, 1]]),
        ("SP2",  [7, 1, "C(O)(=O)(C)",  0, "", ["O", 0, 1]]),
        ("P2",   [7, 1, "C(O)(=O)(C=C)",0, "", ["O", 0, 1]]),
        # ester
        ("N4a",  [7, 0, "COC=O",        1, "", []]),
        ("SN4a", [7, 0, "OC=O",         1, "", ["O", 0, 0]]),
        ("N4a",  [7, 1, "C(O)(=O)(CC)", 0, "", ["O", 0, 0]]),
        ("SN4a", [7, 1, "C(O)(=O)(C)",  0, "", ["O", 0, 0]]),
        ("N4a",  [7, 1, "C(O)(=O)(C=C)",0, "", ["O", 0, 0]]),
        ("N4a",  [7, 1, "C(OC)(=O)(C)", 0, "", []]),
        # Acetal/ketal
        ("SN4a", [7, 0, "OCO",          0, "", ["O", 0, 0]]),
        ("N4a",  [7, 0, "OCOC",         0, "", ["O", 0, 0]]),
        ("SN4a", [7, 1, "C(O)(O)(C)",   0, "", ["O", 0, 0]]), 
        ("N4a",  [7, 1, "C(O)(O)(CC)",  0, "", ["O", 0, 0]]), 
        ("N4a",  [7, 1, "C(O)(OC)(C)",  0, "", ["O", 0, 0]]),
        ("N4a",  [7, 1, "C(O)(OC)(=C)", 0, "", ["O", 0, 0]]),
        # hemiacetal/hemiketal
        ("SP2",  [7, 0, "OCO",          0, "", ["O", 0, 1]]),
        ("P2",   [7, 0, "OCOC",         0, "", ["O", 0, 1]]),
        ("SP2",  [7, 1, "C(O)(O)(C)",   0, "", ["O", 0, 1]]), 
        ("P2",   [7, 1, "C(O)(O)(CC)",  0, "", ["O", 0, 1]]), 
        ("P2",   [7, 1, "C(O)(OC)(C)",  0, "", ["O", 0, 1]]),
        ("P2",   [7, 1, "C(O)(OC)(=C)", 0, "", ["O", 0, 1]]),
        # diol
        ("P4",   [7, 0, "OCO",          0, "", ["O", 1, 1]]),
        ("P4",   [7, 1, "C(O)(O)(C)",   0, "", ["O", 1, 1]]), 
        ("P4",   [7, 1, "C(O)(O)(CC)",  0, "", ["O", 1, 1]]), 
        ("P4",   [7, 0, "OCCO",         0, "", ["O", 1, 1]]), 
        # primary amide
        ("SP5",  [7, 0, "NC=O",         1, "", ["N", 2]]), #SP3a
        ("SP5",  [7, 1, "C(N)(=O)(C)",  0, "", ["N", 2]]), #SN5
        ("P5",   [7, 1, "C(N)(=O)(CC)", 0, "", []]),
        ("P5",   [7, 1, "C(N)(=O)(C=C)",0, "", []]),
        # secondary amide
        ("SP3",  [7, 0, "NC=O",         1, "", ["N", 1]]),
        ("SP3",  [7, 1, "C(N)(=O)(C)",  0, "", ["N", 1]]),
        ("P3",   [7, 0, "O=CNC",        0, "", ["N", 1]]), #P1
        ("P3",   [7, 1, "C(NC)(=O)(C)", 0, "", ["N", 1]]), #P4
        # tertiary amide
        ("SP3a", [7, 0, "NC=O",         1, "", ["N", 0]]),
        ("SP3a", [7, 1, "C(N)(=O)(C)",  0, "", ["N", 0]]), 
        ("P3a",  [7, 0, "O=CNC",        0, "", ["N", 0]]), 
        ("P3a",  [7, 1, "C(NC)(=O)(C)", 0, "", ["N", 0]]),
        ("P3a",  [7, 1, "N(C=O)(C)(C)", 0, "", []]), #P4a
        # sulfoxide
        ("P4",   [7, 1, "S(=O)(C)(C)",    0, "", []]),
        ("SP4",  [7, 1, "CS=O",           0, "", []]), #SP6
        # amino acid
        ("P6",   [7, 1, "C(CN)(=O)(O)",   0, "", []]),
        # primary organoiodine
        ("X1",  [7, 0, "CCI",          2, "", []]),
        ("X1",  [7, 0, "CCCI",         2, "", []]),
        ("SX1", [7, 0, "CI",           2, "", []]),
        # primary organobromide
        ("X2",  [7, 0, "CCBr",         2, "", []]), #SX1
        ("X2",  [7, 0, "BrCC",         2, "", []]), #SX1
        ("X2",  [7, 0, "CCCBr",        2, "", []]), #new
        ("X2",  [7, 0, "BrCCC",        2, "", []]), #new
        ("SX2", [7, 0, "CBr",          2, "", []]),
        ("SX2", [7, 0, "BrC",          2, "", []]),
        # primary organochloride
        ("SX3", [7, 0, "CCCl",         2, "", []]),
        ("SX3", [7, 0, "ClCC",         2, "", []]),
        ("X3",  [7, 0, "CCCCl",        2, "", []]), #new
        ("X3",  [7, 0, "ClCCC",        2, "", []]), #new
        ("TX3", [7, 0, "CCl",          0, "", []]),
        ("TX3", [7, 0, "ClC",          0, "", []]),
        ("X3 ", [7, 1, "C(C)(C)(CCl)", 0, "", []]),
        # primary organoflourine
        ("X4e",  [7, 0, "CCCF",         2, "", []]), #new
        ("SX4e", [7, 0, "CCF",          2, "", []]),
        ("TX4e", [7, 0, "CF",           2, "", []]), #TX3e
        # dichloroethane
        ("SX3h", [7, 1, "C(C)(Cl)(Cl)",0, "", []]),
        ("SX3h", [7, 0, "ClCCl",        0, "", []]),
        ("X3h",  [7, 0, "ClCCCl",       0, "", []]),
        # tetrachloromethane
        ("X1",   [7, 2, "C(Cl)(Cl)(Cl)(Cl)", 0, "."]), #X2
        # sulfonate
        ("SQ4n", [7, 1, "S(=O)(=O)(O)",   0, "", []]), #SP6
        
        # INFERRED AND TUNED from old model
        # di-halogen-o-ethane
        ("SX4h", [7, 0, "FCF",          1, "", []]),
        ("X4h",  [7, 0, "FCCF",         1, "", []]),
        ("X2h",  [7, 0, "BrCCBr",       0, "", []]),
        ("SX2h", [7, 0, "BrCBr",        0, "", []]),
        # glycolaldehyde
        ("N5a",  [7, 1, "C(C)(=O)(CO)", 0, "", []]),
        ("N5a",  [7, 0, "O=CCO",        0, "", []]),
        # Carbamate
        ("SP2",  [7, 1, "C(O)(=O)(N)",  0, "", []]),
        ("P6",   [7, 1, "C(NC)(=O)(O)", 0, "", []]),
        # hydrazone
        ("N1a",  [7, 1, "C(C)(C)(=NN)", 0, "", []]),
        ("N1a",  [7, 0, "CC=NN",        0, "", []]),
        ("SN1a", [7, 0, "C=NN",         0, "", []]),
        # secondary organo-halogens
        ("SX4e", [7, 1, "C(C)(C)(F)",      0, "", []]),
        ("X1",   [7, 1, "C(C)(C)(I)",      0, "", []]),
        ("SX2",  [7, 1, "C(C)(C)(Br)",     0, "", []]),
        ("SX3",  [7, 1, "C(C)(C)(Cl)",     0, "", []]),
        # tertiary organo-halogens
        ("X3",   [7, 2, "C(C)(C)(C)(Cl)"   , 0, "", []]),
        ("X2",   [7, 2, "C(C)(C)(C)(Br)"   , 0, "", []]),
        ("X1",   [7, 2, "C(C)(C)(C)(I)"    , 0, "", []]),
        # tetra-halogen-methane
        ("X2",   [7, 2, "C(Br)(Br)(Br)(Br)", 0, "."]), #X1
        ("X3",   [7, 2, "C(F)(F)(F)(F)"    , 0, "."]), #X4

        # -------------------------------------------
        # Section 8: node that is in 2 benzene rings
        # -------------------------------------------
        
        # FROM THE BUILDING BLOCK TABLE: (unchangable)
        ("TC5e", [8, 2, "CC", 4, "", []]),
        
        # FROM TABLE 24: (no change)
        # aromatics
        ("TC5",  [8, 1, "CC", 4, "", []]),

        # --------------------------------
        # Section 9: 3/4/7/8‑rings section
        # --------------------------------
        # FROM TABLE 24: (no change)
        # cyclic alkane
        ("SC3",  [9, 0, "CCC", 2, "", []]),
        # diol
        ("P4",   [9, 0, "COC", 2, "", []]), #SN4a
        
        # ----------------------------------
        # Section 10: Singular atom section
        # ----------------------------------
        # INFERRED AND TUNED from old model
        # Thiophene
        ("TC6",  [10, 0, "S", 0, "", []]),
        
        # ---------------------------------------
        # Section 11: Structures that are forced
        # ---------------------------------------
        
        # INFERRED and yet to be tuned
        ("P3?",   [11, 2, "S(C)(=O)(=O)(C)"  , 0, "", []]),
        # Tetramethyl-group 14
        ("C2?",   [11, 2, "Ge(C)(C)(C)(C)"   , 0, "", []]),
        ("C2?",   [11, 2, "Pb(C)(C)(C)(C)"   , 0, "", []]),
        ("C2?",   [11, 2, "Si(C)(C)(C)(C)"   , 0, "", []]),
        ("C2?",   [11, 2, "Sn(C)(C)(C)(C)"   , 0, "", []]),
        # sulfate
        ("SQ4n?", [11, 2, "S(=O)(=O)(O)(O)",   0, "", []]),
        ("SC2?",  [11, 1, "Si(C)(C)(C)",    0, "", []]),
        ("P2?",   [11, 1, "C(O)(=O)(CO)",   0, "", []]),
        ("SX4e?", [11, 1, "C(Cl)(F)(F)",    1, "", []]),
        ("X4e?",  [11, 1, "C(OC)(F)(F)",    1, "", []]),
        ("X4e?",  [11, 1, "C(CO)(F)(F)",    1, "", []]),
        ("SX4e?", [11, 1, "C(O)(F)(F)",     1, "", []]),
        ("X4e?",  [11, 1, "C(=C)(F)(F)",    1, "", []]),
        ("X4e?",  [11, 1, "C(=CF)(F)(F)",   1, "", []]),
        ("X4e?",  [11, 1, "C(Cl)(CF)(F)",   1, "", []]),
        ("X4e?",  [11, 1, "C(=O)(CF)(N)",   1, "", []]),
        ("X2?",   [11, 1, "C(F)(Cl)(Cl)",   0, "", []]),
        ("X2?",   [11, 1, "C(CCl)(Cl)(Cl)", 0, "", []]),
        ("SX2?",  [11, 1, "C(C)(Cl)(Cl)",   0, "", []]),
        ("SX2?",  [11, 1, "C(C)(Cl)(CCl)",  0, "", []]),
        ("SX2?",  [11, 1, "C(=C)(Cl)(Cl)",  0, "", []]),
        ("X1?",   [11, 1, "C(Br)(Cl)(Br)",  0, "", []]),
        ("X1?",   [11, 1, "C(Br)(Br)(Br)",  0, "", []]),
        ("X1?",   [11, 1, "C(=CBr)(Br)(Br)",0, "", []]),
        ("X2?",   [11, 1, "C(Cl)(Cl)(Br)",  0, "", []]),
        ("SP2?",  [11, 1, "C(C)(=O)(=O)",   0, "", []]),
        ("N5a?",  [11, 1, "C(C)(C=O)(=C)",  0, "", []]),
        ("N5?",   [11, 1, "C(=C)(C)(CO)",   0, "", []]),
        ("P4?",   [11, 1, "C(O)(C)(CO)",    0, "", []]),
        ("SN4?",  [11, 1, "N(O)(C)(C)",     0, "", []]),
        ("SN4?",  [11, 1, "N(N=O)(C)(C)",   0, "", []]),
        ("SP4?",  [11, 1, "C(O)(=O)(=O)",   0, "", []]),
        ("P5?",   [11, 1, "C(CN)(=O)(N)",   0, "", []]),
        ("P6?",   [11, 1, "C(OC)(=O)(N)",   0, "", []]),
        # Urea
        ("P5?",   [11, 1, "C(NC)(=O)(N)",   0, "", []]),
        ("SP5?",  [11, 1, "C(N)(=O)(N)",    0, "", []]),
        ("P6?",   [11, 1, "C(NO)(=O)(N)",   0, "", []]),
        ("P5?",   [11, 1, "C(NO)(=O)(C)",   0, "", []]),
        ("P5?",   [11, 1, "C(NN)(=O)(C)",   0, "", []]),
        ("SP5?",  [11, 1, "C(N)(=S)(N)",    0, "", []]),
        ("P6?",   [11, 1, "C(N)(=O)(NN)",   0, "", []]),
        ("SP5?",  [11, 1, "C(N)(=N)(N)",    0, "", []]),
        ("X3?",   [11, 1, "C(O)(C)(CCl)",   0, "", []]),
        ("SX3?",  [11, 1, "C(C)(=C)(Cl)",   0, "", []]),
        ("X3?",   [11, 1, "C(C)(CC)(Cl)",   0, "", []]),
        ("X2?",   [11, 1, "C(C)(C)(CBr)",   0, "", []]),
        ("SC3?",  [11, 1, "Pb(C)(C)(C)",    0, "", []]),
        ("SN4a?", [11, 0, "NCO",           1, "", []]),
        ("TP4?",  [11, 0, "OO",            0, "", []]),
        ("TN3a?", [11, 0, "NO",            0, "", []]),
        ("TN3a?", [11, 0, "N=O",           0, "", []]),
        ("TN1a?", [11, 0, "NN",            0, "", []]),
        ("TC5?",  [11, 0, "N=N",           0, "", []]),
        ("SN2?",  [11, 0, "NC=S",          0, "", []]),
        ("SN2?",  [11, 0, "N=C=S",         0, "", []]),
        ("SN2?",  [11, 0, "NC=N",          0, "", []]),
        ("SP1?",  [11, 0, "C=C=O",         0, "", []]),
        ("TC6?",  [11, 0, "C=S",           0, "", []]),
        ("SC6?",  [11, 0, "CS=O",          0, "", []]),
        ("SC6?",  [11, 0, "S=C=S",         0, "", []]),
        ("SN2?",  [11, 0, "S=PC",          0, "", []]),
        ("SX2?",  [11, 0, "C=CBr",         2, "", []]),
        ("SX2?",  [11, 0, "BrCBr",         0, "", []]),
        ("SX3h?", [11, 0, "ClCF",          0, "", []]),
        ("SX3h?", [11, 0, "ClHg",          0, "", []]),
        ("SX3h?", [11, 0, "HgCl",          0, "", []]),
        ("SX3h?", [11, 0, "ClN",           0, "", []]),
        ("SX3h?", [11, 0, "NCl",           0, "", []]),
        ("SX3?",  [11, 0, "ClCO",          0, "", []]),
        ("SC3?",  [11, 0, "CPbC",          0, "", []]),
        ("SC3?",  [11, 0, "CSiC",          0, "", []]),
    ]

    # 2) Create the final dictionary, inserting each entry in order.
    #    If a key already exists, append “+” until it’s unique.
    martini_dict = {}
    for key, val in entries:
        new_key = key
        while new_key in martini_dict:
            new_key = new_key + "+"
        # Use a copy of val (so mutating one list doesn’t affect others, if ever)
        martini_dict[new_key] = list(val)

    return martini_dict
