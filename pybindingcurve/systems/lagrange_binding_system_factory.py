from inspect import signature
from copy import deepcopy

class LagrangeBindingSystemFactory():
    """
    LagrangeBindingSystemFactory produces custom lagrange binding fucntions

    From a simple definition string, such as:
    p+l<->pl
    p+i<->pi

    Parameters
    ----------
    bindingsystem : func
        A function deriving the concentration of complex formed. This fuction 
        may be varied based on different protein-ligand binding systems.
    analytical: bool
        Perform a analytical analysis (default = False)
    """
    default_readout = None
    species = []
    fundamental_species=[]
    reactions_dictionary=None
    original_reactions_dictionary=None
    fundamental_species_in_products=None
    func_string=None
    _add_nonzero_constraints=False
    
    def __init__(self, system_string, add_nonzero_constraints=False):
        
        # Adding non-zero constraints used in testing and development
        self._add_nonzero_constraints=add_nonzero_constraints
        
        # Get reactions tuples
        reactions = self._get_reactions_and_set_readout(system_string.lower())
        # Populate self.species, an ordered list of species encountered
        [self.species.append(s) for r in reactions for s in r if not s in self.species]
        
        # Make reactions dictionary and then simplify it
        self.original_reactions_dictionary={k[2]:[] for k in reactions}
        [self.original_reactions_dictionary[p].append([r1,r2,self._get_correct_kd((r1,r2,p))]) for r1,r2,p in reactions]

        self.reactions_dictionary=self._simplify_reactions_dictionary(self.original_reactions_dictionary)
        
        # Find fundamental species, ordered by occurrence in systems definition string
        self.fundamental_species=[x for x in self.species if x in set(self.species)-set(self.reactions_dictionary.keys())]
        
        assert len(
            self.fundamental_species) > 0, "Malformed system, no fundamental species"
        assert len(self.reactions_dictionary.keys()) > 0, "Malformed species, no products"
        # Find fundamental species counts in products
        self.fundamental_species_in_products={}
        for product,counts in self._get_num_fundamental_species_by_products(self.original_reactions_dictionary).items():
            for fundamental_species,count in counts.items():
                if fundamental_species not in self.fundamental_species_in_products.keys():
                    self.fundamental_species_in_products[fundamental_species]={}
                if product not in self.fundamental_species_in_products[fundamental_species].keys():
                    self.fundamental_species_in_products[fundamental_species][product]=count
                else:
                    self.fundamental_species_in_products[fundamental_species][product]+=count

        self.func_string=self.get_func_string()
        exec(self.func_string, globals())
        self.binding_function=eval("custom_lagrange_binding_system")
        self.custom_function_arguments=list(signature(self.binding_function).parameters.keys())

    def write_func_to_python_file(self, filename):
        out_file=open(filename, "w")
        out_file.write(self.func_string)
        out_file.close()

    def get_func_string(self):
        # Write header and function definition
        custom_lagrange_definition = '\"\"\"\nCustom generated binding system\n\nLagrane multiplier binding system genetated with generated with \nhttps://github.com/stevenshave/lagrange-binding-systems/write_custom_system.py\"\"\"\n\n'
        custom_lagrange_definition += f'from scipy.optimize import fsolve\nfrom autograd import grad\ndef custom_lagrange_binding_system('
        custom_lagrange_definition += ", ".join([f"{x}0" for x in self.fundamental_species])+", "+", ".join([x[2] for _,r in self.reactions_dictionary.items() for x in r])+"):\n"
        custom_lagrange_definition += "\tdef F(X): # Augmented Lagrange function\n"
        # Write fundamental species concs
        custom_lagrange_definition += "".join([f"\t\t{x}=X[{ix}]\n" for ix, x in enumerate(self.fundamental_species)])
        
        # Write the mass balances
        mass_balances=[]
        for k,v in self.original_reactions_dictionary.items():
            line=f"\t\t{k}=("+"+".join([f"({sr[0]}*{sr[1]})" for sr in v])+")/("
            line+="+".join([sr[2] for sr in v])+")"
            mass_balances.append(line+"\n")
        custom_lagrange_definition += "".join(mass_balances)


        long_mass_balances=[]
        for k,v in self.reactions_dictionary.items():
            line=f"\t\t{k}=("+"+".join([f"({sr[0]}*{sr[1]})" for sr in v])+")/("
            line+="+".join([sr[2] for sr in v])+")"
            long_mass_balances.append(line+"\n")


        # Write constraints
        constraints=[]
        for constraint_num in range(1,len(self.fundamental_species)+1):
            line=f"\t\tconstraint{constraint_num}={self.fundamental_species[constraint_num-1]}0-({self.fundamental_species[constraint_num-1]}+"
            line+="+".join([f"{count}*{product}" for product,count in self.fundamental_species_in_products[self.fundamental_species[constraint_num-1]].items()])+")"
            constraints.append(line+"\n")
        if self._add_nonzero_constraints:
            constraints.append(f"\t\tnonzero_constraint={'-'.join([f'(constraint{x}-abs(constraint{x}))' for x in range(1, len(self.fundamental_species)+1)])}\n")
        custom_lagrange_definition+="".join(constraints)
        
        # Write return statement
        return_statements = ["\t\treturn "+self.default_readout+"-"]
        return_statements.append("+".join([f"X[{i+len(self.fundamental_species)}]*constraint{i+1}" for i in range(len(self.fundamental_species))]))
        if self._add_nonzero_constraints:
            return_statements.append(f"+X[{len(self.fundamental_species)+len(self.fundamental_species)}]*nonzero_constraint")
        return_statements[-1] = return_statements[-1]+"\n"
        custom_lagrange_definition+="".join(return_statements)

        # Write derivation
        custom_lagrange_definition += "\tderivative_function = grad(F) # Gradients of the Lagrange function\n"

        custom_lagrange_definition += "\t"+", ".join([f"{s}" for s in self.fundamental_species])+", "+", ".join([f"lam{x}" for x in range(len(self.fundamental_species)+self._add_nonzero_constraints)]) + \
            " = fsolve(derivative_function, ["+", ".join([f"{x}0" for x in self.fundamental_species]) + \
            ", " + \
            ", ".join(["1.0" for x in range(
                len(self.fundamental_species)+self._add_nonzero_constraints)])+"])\n"
        custom_lagrange_definition += "\treturn {"+", ".join([f"'{s}':{s}" for s in self.fundamental_species])+", "+", ".join(
            [f"'{mb.split('=')[0].strip()}':{mb.split('=')[1].strip()}" for mb in long_mass_balances])+"}\n"

        return custom_lagrange_definition
    def _get_num_fundamental_species_by_products(self, original_reactions_dict):
        fundamental_species_in_products={}
        for product in original_reactions_dict.keys():
            fundamental_species_in_products[product]={}
            for species in self.species:
                c=original_reactions_dict[product][0][0:2].count(species)
                if c>0:
                    fundamental_species_in_products[product][species]=c
        # Counts could contain non-fundamental products at this stage, we need to iteratively simplify.
        while True:
            something_changed=False
            for product, species_dict in fundamental_species_in_products.items():
                non_fundamental=[x for x in species_dict.keys() if x not in self.fundamental_species]
                if len(non_fundamental)==0:continue
                something_changed=True
                for nfs in non_fundamental:
                    del species_dict[nfs]
                    for k, v in fundamental_species_in_products[nfs].items():
                        if k in species_dict.keys():
                            species_dict[k]+=v
                        else:
                            species_dict[k]=v
            if not something_changed:break
        return fundamental_species_in_products

    def _extract_species_numbers(self, string, species):
        loc=-1
        res=[]
        while True:
            loc=string.find(species,loc+1)
            if loc==-1:break
            i=loc+len(species)
            while i<len(string) and string[i].isdigit():
                i+=1
                if i==len(string):break
            res.append(string[loc+len(species):i])
        return res            


    def _get_correct_kd(self, reaction):
        """
        Get the correct KD for a reaction

        Setting the correct KD is complex when we have a system such as:
        P+L<->PL1, P+L<->PL2, PL1+L<->PL1L2, PL2+L<->PL1L2

        """
        r1=reaction[0]
        r2=reaction[1]
        p=reaction[2]

        #If no digits, just return the KD term
        if (not(any(map(lambda x:x.isdigit(), r1)))) and (not(any(map(lambda x:x.isdigit(), r2)))) and (not(any(map(lambda x:x.isdigit(), p)))):
            return f"kd_{r1}_{r2}"
        
        # There are digits present.
        species_numbers_r2_in_r1=self._extract_species_numbers(r1, r2)
        species_numbers_r1_in_r2=self._extract_species_numbers(r2, r1)
        species_numbers_r1_in_p=self._extract_species_numbers(p, r1)
        species_numbers_r2_in_p=self._extract_species_numbers(p, r2)

        missing_from_r1=list(set(species_numbers_r2_in_p)-set(species_numbers_r2_in_r1))
        missing_from_r2=list(set(species_numbers_r1_in_p)-set(species_numbers_r1_in_r2))

        if '' in missing_from_r1: missing_from_r1.remove('')
        if '' in missing_from_r2: missing_from_r2.remove('')

        assert len(missing_from_r1)+len(missing_from_r2)==1, "A reaction should add one species to another in a clearly understandable way"

        if len(missing_from_r1)>0:
            return f"kd_{r1}_{r2}{next(iter(missing_from_r1))}"
        else:
            return f"kd_{r1}{next(iter(missing_from_r2))}_{r2}"
       

    def _simplify_reactions_dictionary(self, reactions_dict):
        # Iteratively simplify the reactions dictionary until all is expressed
        # with fundamental species
        simplified_reactions_dict=deepcopy(reactions_dict)
        while True:
            something_changed=False
            for k,v in simplified_reactions_dict.items():
                for i_r, r in enumerate(v):
                    for tuple_it in range(2):
                        if r[tuple_it] in simplified_reactions_dict.keys():
                            something_changed=True
                            replacement=simplified_reactions_dict[r[tuple_it]][0]
                            simplified_reactions_dict[k][i_r][tuple_it]=f"{replacement[0]}*{replacement[1]}/{replacement[2]}"
            if not something_changed: break
        return simplified_reactions_dict


    def _get_reactions_and_set_readout(self, system_string):
        """
        Parse a system string to list of reaction tuples

        Gets list of reaction tuples, also sets self.readout if * is found on
        a species.  If no * is found, readout is the product of the first
        reaction.

        """

        # Set up reaction_strings and tuples
        reaction_strings = [item.strip() for sublist in [r.split(
            "\n") for r in system_string.split("#")[0].split(",")] for item in sublist if len(item.strip()) > 0]
        reaction_tuples = [(s[0], s[1].split("<->")[0], s[1].split("<->")[1])
                           for s in [r.split("+") for r in reaction_strings]]

        # Find the readout
        for i,v in enumerate(reaction_tuples):
            print(v)
            if v[0].find("*")>-1:
                self.default_readout=v[0].replace("*","")
                reaction_tuples[i]=(v[0].replace("*",""),v[1], v[2])
            if v[1].find("*")>-1:
                self.default_readout=v[1].replace("*","")
                reaction_tuples[i]=(v[0],v[1].replace("*",""), v[2])
            if v[2].find("*")>-1:
                self.default_readout=v[2].replace("*","")
                reaction_tuples[i]=(v[0],v[1],v[2].replace("*",""))
        if not self.default_readout:
            self.default_readout = reaction_tuples[0][2]

        return reaction_tuples


if __name__ == "__main__":
    custom_system = """
        P+P<->PP*
        P+L<->PL
        PP+L<->PPL1
        PP+L<->PPL2
        PPL1+L<->PPL1L2
        PPL2+L<->PPL1L2
"""
    new_lagrange = LagrangeBindingSystemFactory(custom_system)
    print(f"{new_lagrange.custom_function_arguments}")
    print(new_lagrange.binding_function(10,10,55,10,10,10,10,10))
