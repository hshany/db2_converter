from rdkit import Chem
from db2_converter.utils.utils import next_mol2_lines, ATOMTYPE, BONDTYPE

acceptable_elements = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"]


def infopart(mol2block):
    p1 = mol2block.index("@<TRIPOS>ATOM\n")
    p2 = mol2block.index("@<TRIPOS>BOND\n")
    Natom, Nbond = mol2block[2].strip().split()[0:2]
    startpart = mol2block[:p1]
    atompart = ["@<TRIPOS>ATOM\n"]
    bondpart = ["@<TRIPOS>BOND\n"]
    # atompart (antechamber mol2 format)
    for i in range(p1 + 1, p1 + 1 + int(Natom)):
        items = mol2block[i].strip().split()
        x, y, z = items[2:5]
        chg = items[-1]
        atompart.append(
            ATOMTYPE.format(
                int(items[0]),
                items[1],
                float(x),
                float(y),
                float(z),
                items[5],
                items[6],
                items[7],
                float(chg),
            )
        )
    # bondpart (antechamber mol2 format)
    for i in range(p2 + 1, p2 + 1 + int(Nbond)):
        items = mol2block[i].strip().split()
        bondpart.append(
            BONDTYPE.format(int(items[0]), int(items[1]), int(items[2]), items[3])
        )
    return startpart, atompart, bondpart


def fixoption(mol, mol2block, options):
    Nitro = Chem.MolFromSmarts("[N+](=O)[O-]")  # NO2-
    Phos1 = Chem.MolFromSmarts("P(=O)([O-])[O-]")  # PO3 2-
    Phos2 = Chem.MolFromSmarts("P(=O)([O-])")  # PO2-
    COO = Chem.MolFromSmarts("O=C[O-]")  # COO-
    dithio = Chem.MolFromSmarts("C(=S)-S")
    oxidon = Chem.MolFromSmarts("[O-][n+]")

    if options.fixnitroante:  # after antechamber
        tgt_moltype = "N.pl3"
        groups = [Nitro]
        mol2block = fix_nitroante(mol, mol2block, tgt_moltype, groups)
    if options.fixnitro:
        tgt_moltype = "N.2"
        groups = [Nitro]  # NO2-
        mol2block = fix_nitro_phos_coo(mol, mol2block, tgt_moltype, groups)
    if options.fixphos:
        tgt_moltype = "P.3"
        groups = [Phos2, Phos1]  # PO2 -, PO3 2- ! change order to ensure largest match
        mol2block = fix_nitro_phos_coo(mol, mol2block, tgt_moltype, groups)
    if options.fixcoo:
        tgt_moltype = "C.2"
        groups = [COO]  # COO-
        mol2block = fix_nitro_phos_coo(mol, mol2block, tgt_moltype, groups)
    if options.fixdithio:  # after antechamber
        groups = [dithio]
        mol2block = fix_dithioic(mol, mol2block, groups)
    if options.fixoxido:
        groups = [oxidon]
        mol2block = fix_oxidopyridine(mol, mol2block, groups)

    return mol2block


def fix_nitroante(mol, mol2block, tgt_moltype, groups):
    startpart, atompart, bondpart = infopart(mol2block)
    for group in groups:
        matches = mol.GetSubstructMatches(group)
    for match in matches:
        Noxys = []
        Pidx = -1
        for idx in match:
            # idx starts from 0, correspondingly atompart 1
            chg = atompart[idx + 1].strip().split()[-1]
            mol2type = atompart[idx + 1].strip().split()[5]
            if mol2type == tgt_moltype:
                Pidx = idx
            if mol2type == "O.2":
                Noxys.append(idx)
        if Pidx == -1:
            break  # not match
        for i, bond in enumerate(bondpart[1:]):
            if bond != "\n":
                bondinfo = bond.strip().split()
                bondorder = bondinfo[3]
                if int(bondinfo[1]) == Pidx + 1 or int(bondinfo[2]) == Pidx + 1:
                    bondorder = str(bondorder).ljust(4)
                    if int(bondinfo[1]) - 1 in Noxys or int(bondinfo[2]) - 1 in Noxys:
                        bondorder = "2   "
                    items = bondpart[i + 1].strip().split()
                    bondpart[i + 1] = BONDTYPE.format(
                        int(items[0]), int(items[1]), int(items[2]), bondorder
                    )

    return startpart + atompart + bondpart + ["\n"]


def fix_nitro_phos_coo(mol, mol2block, tgt_moltype, groups):
    startpart, atompart, bondpart = infopart(mol2block)
    for group in groups:
        matches = mol.GetSubstructMatches(group)
        for match in matches:
            Pidx = -1
            sb, db = [], []
            for idx in match:
                # idx starts from 0, correspondingly atompart 1
                chg = atompart[idx + 1].strip().split()[-1]
                mol2type = atompart[idx + 1].strip().split()[5]
                if mol2type == tgt_moltype:
                    Pidx = idx
                if mol2type == "O.co2":
                    if int(float(chg)) == -1:
                        sb.append(idx)
                    elif int(float(chg)) == 0:
                        db.append(idx)
            for i, bond in enumerate(bondpart[1:]):
                if bond != "\n":
                    bondinfo = bond.strip().split()
                    bondorder = bondinfo[3]
                    BOindex = bond.rindex(bondorder)
                    if int(bondinfo[1]) == Pidx + 1 or int(bondinfo[2]) == Pidx + 1:
                        bondorder = str(bondorder).ljust(4)
                        if int(bondinfo[1]) - 1 in sb or int(bondinfo[2]) - 1 in sb:
                            atomidx = (
                                int(bondinfo[1])
                                if int(bondinfo[2]) == Pidx + 1
                                else int(bondinfo[2])
                            )
                            mol2type = "O.3  "
                            bondorder = "1   "
                        elif int(bondinfo[1]) - 1 in db or int(bondinfo[2]) - 1 in db:
                            atomidx = (
                                int(bondinfo[1])
                                if int(bondinfo[2]) == Pidx + 1
                                else int(bondinfo[2])
                            )
                            mol2type = "O.2  "
                            bondorder = "2   "
                        else:
                            continue
                        if not tgt_moltype == "C.2":  # not change COO-
                            items = atompart[atomidx].strip().split()
                            x, y, z = items[2:5]
                            chg = items[-1]
                            atompart[atomidx] = ATOMTYPE.format(
                                int(items[0]),
                                items[1],
                                float(x),
                                float(y),
                                float(z),
                                mol2type,
                                items[6],
                                items[7],
                                float(chg),
                            )
                        items = bondpart[i + 1].strip().split()
                        bondpart[i + 1] = BONDTYPE.format(
                            int(items[0]), int(items[1]), int(items[2]), bondorder
                        )
    return startpart + atompart + bondpart + ["\n"]


def fix_dithioic(mol, mol2block, groups):
    Cmol2type = "C.2  "
    S1mol2type = "S.2  "  # C=S
    S2mol2type = "S.3  "  # C-S-
    CS1bondorder = "2   "
    CS2bondorder = "1   "
    startpart, atompart, bondpart = infopart(mol2block)
    for group in groups:
        matches = mol.GetSubstructMatches(group)
    for match in matches:
        Cidx, S1idx, S2idx = match
        S1chg = atompart[S1idx + 1].strip().split()[8]
        S2chg = atompart[S2idx + 1].strip().split()[8]
        if int(float(S1chg)) == -1 and int(float(S2chg)) == 0:
            S1idx, S2idx = S2idx, S1idx
        elif int(float(S1chg)) == 0 and int(float(S2chg)) == -1:
            pass
        else:
            print(">>> something unexpected happens!")
            break
        for atomidx in range(len(atompart)):
            if atomidx - 1 == Cidx:
                items = atompart[atomidx].strip().split()
                atompart[atomidx] = ATOMTYPE.format(
                    int(items[0]),
                    items[1],
                    float(items[2]),
                    float(items[3]),
                    float(items[4]),
                    Cmol2type,
                    items[6],
                    items[7],
                    float(items[-1]),
                )
            if atomidx - 1 == S1idx:
                items = atompart[atomidx].strip().split()
                atompart[atomidx] = ATOMTYPE.format(
                    int(items[0]),
                    items[1],
                    float(items[2]),
                    float(items[3]),
                    float(items[4]),
                    S1mol2type,
                    items[6],
                    items[7],
                    float(items[-1]),
                )
            if atomidx - 1 == S2idx:
                items = atompart[atomidx].strip().split()
                atompart[atomidx] = ATOMTYPE.format(
                    int(items[0]),
                    items[1],
                    float(items[2]),
                    float(items[3]),
                    float(items[4]),
                    S2mol2type,
                    items[6],
                    items[7],
                    float(items[-1]),
                )
        for i, bond in enumerate(bondpart[1:]):
            if bond != "\n":
                bondinfo = bond.strip().split()
                bondorder = bondinfo[3]
                if int(bondinfo[1]) == Cidx + 1 or int(bondinfo[2]) == Cidx + 1:
                    bondorder = str(bondorder).ljust(4)
                    if int(bondinfo[1]) == S1idx + 1 or int(bondinfo[2]) == S1idx + 1:
                        bondorder = CS1bondorder
                    if int(bondinfo[1]) == S2idx + 1 or int(bondinfo[2]) == S2idx + 1:
                        bondorder = CS2bondorder
                    # bondpart[i+1] = bondpart[i+1][:19] + bondorder + "\n"
                    items = bondpart[i + 1].strip().split()
                    bondpart.append(
                        BONDTYPE.format(
                            int(items[0]), int(items[1]), int(items[2]), bondorder
                        )
                    )

    return startpart + atompart + bondpart + ["\n"]


def fix_oxidopyridine(mol, mol2block, groups):
    Omol2type = "O.3  "
    Nmol2type = "N.pl3"
    ONbondorder = "1   "
    startpart, atompart, bondpart = infopart(mol2block)
    for group in groups:
        matches = mol.GetSubstructMatches(group)
    for match in matches:
        Oidx, Nidx = match
        for atomidx in range(len(atompart)):
            if atomidx - 1 == Oidx:
                # atompart[atomidx] = atompart[atomidx][:50] + Omol2type + atompart[atomidx][55:]
                items = atompart[atomidx].strip().split()
                atompart[atomidx] = ATOMTYPE.format(
                    int(items[0]),
                    items[1],
                    float(items[2]),
                    float(items[3]),
                    float(items[4]),
                    Omol2type,
                    items[6],
                    items[7],
                    float(items[-1]),
                )
            if atomidx - 1 == Nidx:
                # atompart[atomidx] = atompart[atomidx][:50] + Nmol2type + atompart[atomidx][55:]
                items = atompart[atomidx].strip().split()
                atompart[atomidx] = ATOMTYPE.format(
                    int(items[0]),
                    items[1],
                    float(items[2]),
                    float(items[3]),
                    float(items[4]),
                    Nmol2type,
                    items[6],
                    items[7],
                    float(items[-1]),
                )
        for i, bond in enumerate(bondpart[1:]):
            if bond != "\n":
                bondinfo = bond.strip().split()
                bondorder = bondinfo[3].ljust(4)
                if (int(bondinfo[1]) == Oidx + 1 and int(bondinfo[2]) == Nidx + 1) or (
                    int(bondinfo[2]) == Oidx + 1 and int(bondinfo[1]) == Nidx + 1
                ):
                    bondorder = ONbondorder
                    # bondpart[i+1] = bondpart[i+1][:19] + bondorder + "\n"
                    items = bondpart[i + 1].strip().split()
                    bondpart.append(
                        BONDTYPE.format(
                            int(items[0]), int(items[1]), int(items[2]), bondorder
                        )
                    )

    return startpart + atompart + bondpart + ["\n"]


class FixOptions:
    def __init__(self):
        self.fixnitroante = True
        self.fixphos = True
        self.fixcoo = True
        self.fixnitro = True
        self.fixdithio = False
        self.fixoxido = False


def fixmol2(insmi, inpmol2, options):
    mol2block_fix = []
    mol2file = inpmol2
    all_blocks = [x for x in next_mol2_lines(mol2file)]
    mol = Chem.MolFromSmiles(insmi)
    for mol2block in all_blocks:
        mol2block = fixoption(mol, mol2block, options)
        mol2block_fix.append("".join(mol2block))
    ### overwrite raw .mol2 file
    with open(mol2file, "w") as f:
        f.write("".join(mol2block_fix))

def fixmol2_using_sdf_rdmol(sdfrdmol, inpmol2, options):
    mol2block_fix = []
    mol = sdfrdmol
    mol2file = inpmol2
    all_blocks = [x for x in next_mol2_lines(mol2file)]
    for mol2block in all_blocks:
        mol2block = fixoption(mol, mol2block, options)
        mol2block_fix.append("".join(mol2block))
    ### overwrite raw .mol2 file
    with open(mol2file, "w") as f:
        f.write("".join(mol2block_fix))

def fixmol2_by_template(inpmol2, tempmol2):
    mol2block_fix = []
    tempfile = tempmol2
    tempblock = [x for x in next_mol2_lines(tempfile)][0]
    mol2file = inpmol2
    mol2blocks = [x for x in next_mol2_lines(mol2file)]  # Read ALL conformers
    _, Tatompart, Tbondpart = infopart(tempblock)

    # Process each conformer
    for mol2block in mol2blocks:
        startpart, atompart, bondpart = infopart(mol2block)
        newatompart = []
        for i in range(len(atompart)):
            if (
                i == 0 or len(atompart[i]) < 70
            ):  # @<TRIPOS>ATOM and other potential parts, 70 is arbitrary set but could be suitable for most cases
                newatompart.append(atompart[i])
            else:
                linesplit = atompart[i].strip().split()
                items = Tatompart[i].strip().split()
                x, y, z = linesplit[2:5]
                chg = linesplit[-1]
                newatompart.append(
                    ATOMTYPE.format(
                        int(items[0]),
                        items[1],
                        float(x),
                        float(y),
                        float(z),
                        items[5],
                        items[6],
                        items[7],
                        float(chg),
                    )
                )
        mol2block_fix.append("".join(startpart + newatompart + Tbondpart + ["\n"]))

    ### overwrite raw .mol2 file
    with open(mol2file, "w") as f:
        f.write("".join(mol2block_fix))


def fixmol2_so2_by_template(inpmol2, tempmol2):
    tempmol = Chem.MolFromMol2File(tempmol2, sanitize=False, removeHs=False)
    mol = Chem.MolFromMol2File(inpmol2, sanitize=False, removeHs=False)

    tempmol2lines = list(next_mol2_lines(tempmol2))[0]
    tempATOMlineidx = tempmol2lines.index("@<TRIPOS>ATOM\n")
    tempBONDlineidx = tempmol2lines.index("@<TRIPOS>BOND\n")

    allmol2lines = list(next_mol2_lines(inpmol2))

    so2groups = tempmol.GetSubstructMatches(Chem.MolFromSmarts("S(=O)(=O)[O-]"))

    if so2groups:
        print(">>> so2 groups detected!")
        print(so2groups)
        newmol2lines = []
        for mol2lines in allmol2lines:
            ATOMlineidx = mol2lines.index("@<TRIPOS>ATOM\n")
            BONDlineidx = mol2lines.index("@<TRIPOS>BOND\n")
            for so2group in so2groups:
                for idx in so2group:
                    items = tempmol2lines[idx + tempATOMlineidx + 1].split()  # ATOM
                    olditems = mol2lines[idx + ATOMlineidx + 1].split()  # ATOM
                    x = olditems[2]
                    y = olditems[3]
                    z = olditems[4]
                    chg = olditems[-1]
                    mol2lines[
                        idx + ATOMlineidx + 1
                    ] = f"    {int(items[0]):>3d} {items[1]:<8s} {float(x):>10.4f} {float(y):>10.4f} {float(z):>10.4f} {items[5]:<6s}{items[6]:>6s} {items[7]:<6s}{float(chg):>12.6f}\n"  # according to antechamber mol2 format
                    if mol.GetAtomWithIdx(idx).GetAtomicNum() == 16:
                        Sidx = idx
                for idx in so2group:
                    if idx != Sidx:
                        bondidx = mol.GetBondBetweenAtoms(Sidx, idx).GetIdx()
                        tempbondidx = tempmol.GetBondBetweenAtoms(Sidx, idx).GetIdx()
                        items = tempmol2lines[tempbondidx + tempBONDlineidx + 1].split()
                        olditems = mol2lines[bondidx + BONDlineidx + 1].split()
                        mol2lines[
                            bondidx + BONDlineidx + 1
                        ] = f"{int(olditems[0]):>6d}{int(items[1]):>6d}{int(items[2]):>6d} {items[3]:<4s}\n"  # according to antechamber mol2 format
            newmol2lines.append("".join(mol2lines))
        with open(inpmol2, "w") as f:
            f.write("".join(newmol2lines))


def mol22smi(mol2file):
    try:
        mol = Chem.MolFromMol2File(mol2file)
        smi = Chem.MolToSmiles(mol, isomericSmiles=False)
        return smi
    except:
        return


def fixDuatom(mol2file):
    newmol2blocks = []
    for mol2block in next_mol2_lines(mol2file):
        startpart, atompart, bondpart = infopart(mol2block)
        newatompart = [atompart[0]]
        for line in atompart[1:]:
            items = line.strip().split()
            element = items[5]
            if element.split(".")[0] not in acceptable_elements:
                element = items[1].split(items[0])[0]
            x, y, z = items[2:5]
            chg = items[-1]
            newatompart.append(
                ATOMTYPE.format(
                    int(items[0]),
                    items[1],
                    float(x),
                    float(y),
                    float(z),
                    element,
                    items[6],
                    items[7],
                    float(chg),
                )
            )

        newmol2blocks.append("".join(startpart + newatompart + bondpart + ["\n"]))
    with open(mol2file, "w") as f:
        f.write("".join(newmol2blocks))


def fixmol2_and_du(smi, tmp0fixmol2):
    tmpmol = Chem.MolFromSmiles(smi)
    fixDu = False
    for atom in tmpmol.GetAtoms():
        if atom.GetSymbol not in acceptable_elements:
            fixDu = True
            break
    fixoptions = FixOptions()
    fixmol2(smi, tmp0fixmol2, fixoptions)
    if fixDu:
        fixDuatom(tmp0fixmol2)

def unify_index_order_mol2file(inmol2file,outmol2file):
    # This is beneficial if different mol2blocks in inmol2file does not share the same index order
    all_mol2blocks = []
    for mol2lines in next_mol2_lines(inmol2file):
        ref_mol2block = mol2lines
        startpart, ref_atompart, ref_bondpart = infopart(ref_mol2block)
        refmol = Chem.MolFromMol2Block("".join(mol2lines),removeHs=False)
        _ = Chem.MolToSmiles(refmol)
        smi_to_refmol_index = list(map(int,refmol.GetProp("_smilesAtomOutputOrder")[1:-2].split(",")))
        break
    for mol2lines in next_mol2_lines(inmol2file):
        probe_mol2block = mol2lines
        startpart, probe_atompart, probe_bondpart = infopart(probe_mol2block)
        probemol = Chem.MolFromMol2Block("".join(probe_mol2block),removeHs=False)
        _ = Chem.MolToSmiles(probemol)
        smi_to_probemol_index = list(map(int,probemol.GetProp("_smilesAtomOutputOrder")[1:-2].split(",")))
        probe_ref_atom_dict = dict(zip(smi_to_refmol_index,smi_to_probemol_index))
        update_probe_atompart = ["@<TRIPOS>ATOM\n"]
        for i in range(len(probe_ref_atom_dict)):
            update_probe_atompart.append("    {:>3d}".format(i+1) + probe_atompart[probe_ref_atom_dict[i]+1][7:])
        updated_mol2lines = startpart + update_probe_atompart + ref_bondpart + ["\n"]
        all_mol2blocks.append("".join(updated_mol2lines))
    with open(outmol2file,"w") as f:
        f.write("".join(all_mol2blocks))
