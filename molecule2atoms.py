# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
def parse_molecule (formula):
    
    import re
    #gets the molecules inside the brackets and their bracket valency and also the molecules outside the bracket
    def brckttext_val(txt,brkts):
        start,end=brkts[0],brkts[1]
        txt_outbrkt=txt
        li=[]
        if start in txt:
            pattern=r"(\{}.+?\{})[0-9]*".format(start,end)
            t=re.findall(pattern,txt)
            for each in t:
                start_index=txt_outbrkt.find(each)
                end_index=txt_outbrkt.find(each)+len(each)
                txt_inbrkt=each[1:-1]
                try:    
                    brkt_val=int(txt_outbrkt[end_index])
                    txt_outbrkt=txt_outbrkt[:start_index]+txt_outbrkt[end_index+1:]
                except ValueError:
                    brkt_val=1
                    txt_outbrkt=txt_outbrkt[:start_index]+txt_outbrkt[end_index:]
                li.append((txt_inbrkt,brkt_val))
        else:
            li.append(("",0))
        return txt_outbrkt,li
    #converts invidual molecules to atoms 
    def vals(mols):
        pattern=r"[A-Z][a-z]?[0-9]*"
        mol=re.findall(pattern,mols)
        return mol
    def atom(mol):
        pattern=r"[0-9]+"
        try:
            atoms=re.search(pattern,mol)
            rem=mol[:atoms.start()]+mol[atoms.end():]
            return rem,atoms.group(0)
        except AttributeError:
            return mol,1
    def mols2atoms(mols):
        s=[]
        for each in vals(mols):
            s.append(atom(each))
        return s
    #adds the valency outside the brackets
    def add_outval(lis,outval):
        d=[]
        for each in lis:
            d.append([each[0],int(each[1])*int(outval)])
        return d
    def totalin(b):
        temp=[]
        temp.extend(mols2atoms(b[0]))
        for each in b[1]:
            temp.extend(add_outval(mols2atoms(each[0]),each[1]))
        return temp
    #loop for the complex molecules to change into simpler molecules
    def outinbrackets(molecule,degree=3,li=[]):
        brkt={3:"{}",2:"[]",1:"()"}
        if degree==1:
            li.append(totalin(brckttext_val(molecule,brkt[degree])))
            return li
        (x,[(y,z)])=brckttext_val(molecule,brkt[degree])

        while degree>0:
            out_brkt=outinbrackets(x,degree-1,[])
            in_brkt=outinbrackets(y,degree-1,[])
            li.extend([out_brkt[0]+add_outval(in_brkt[0],int(z))])
            degree=0
    
        return li  
    #changing molecule counter to dictionary
    def list_to_dict(temp):    
        final={}
        for (k,v) in temp:
            if k in final:
                final[k]+=int(v)
            else:
                final[k]=int(v)
        return final
    return list_to_dict(outinbrackets(formula)[0])
if __name__=="__main__":
    #molecule=input("Enter molecule to parse")
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("Molecule_formula", help="display a square of a given number",
                   type=str)
    args = parser.parse_args()
    molecule=args
    print(parse_molecule(args.Molecule_formula))
    
