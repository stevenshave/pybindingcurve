from inspect import signature
import numpy as np
from sys import version_info
from mpmath import mpf, mp, findroot, almosteq

class MinimizerBindingSystemFactory:
	"""
	MinimizerBindingSystemFactory produces custom minimizer-based binding
	system fucntions

	From a simple definition string representing standard 1:1:1 competition, such as in Example 1:

	Example 1:
	----------
	p+l<->pl
	p+i<->pi

	To more complex examples such as example 2:

	Example 2:
	----------
	P+P<->PP*
	P+L<->PL
	PP+L<->PPL1
	PP+L<->PPL2
	PPL1+L<->PPL1L2
	PPL2+L<->PPL1L2

	Example 2 should produce the following function - which is annotated to indicate
	sections which this class builds through parsing of the system definition string.
	Numbers like ### 1 ### denote sections which are described bellow and documented
	in the code

	def custom_minimizer_system(p,l,kd_p_p_pp,kd_p_l_pl,kd_pp_l_ppl1,kd_pp_l_ppl2,kd_ppl1_l_ppl1l2,kd_ppl2_l_ppl1l2): ### 1 ###
			p=mpf(p) ### 2 ###
			l=mpf(l) ### 2 ###
			kd_p_p_pp=mpf(kd_p_p_pp) ### 3 ###
			kd_p_l_pl=mpf(kd_p_l_pl) ### 3 ###
			kd_pp_l_ppl1=mpf(kd_pp_l_ppl1) ### 3 ###
			kd_pp_l_ppl2=mpf(kd_pp_l_ppl2) ### 3 ###
			kd_ppl1_l_ppl1l2=mpf(kd_ppl1_l_ppl1l2) ### 3 ###
			kd_ppl2_l_ppl1l2=mpf(kd_ppl2_l_ppl1l2) ### 3 ###
			def f(p_f,l_f): ### 4 ###
					pp=p_f*p_f/kd_p_p_pp ### 5 ###
					pl=p_f*l_f/kd_p_l_pl ### 5 ###
					ppl1=pp*l_f/kd_pp_l_ppl1 ### 5 ###
					ppl2=pp*l_f/kd_pp_l_ppl2 ### 5 ###
					ppl1l2=(ppl1*l_f+ppl2*l_f)/(kd_ppl1_l_ppl1l2+kd_ppl2_l_ppl1l2) ### 5 ###
					return p0-(p+2*pp+2*ppl1+2*ppl2+pl+2*ppl1l2),l0-(l+ppl1+ppl2+pl+2*ppl1l2) ### 6 ###
			p_f,l_f=findroot(f, [mpf(0), mpf(0)], tol=1e-10) ### 7 ###
			pp=p_f*p_f/kd_p_p_pp ### 5 ###
			pl=p_f*l_f/kd_p_l_pl ### 5 ###
			ppl1=pp*l_f/kd_pp_l_ppl1 ### 5 ###
			ppl2=pp*l_f/kd_pp_l_ppl2 ### 5 ###
			ppl1l2=(ppl1*l_f+ppl2*l_f)/(kd_ppl1_l_ppl1l2+kd_ppl2_l_ppl1l2) ### 5 ###
			return {'p_f':p_f,'l_f':l_f,'pp':pp,'pl':pl,'ppl1':ppl1,'ppl2':ppl2,'ppl1l2':ppl1l2} ### 8 ###

	Section ### 1 ### :	Define the custom minimizer function which takes arguments for fundamental species 
						centration and KDs.
	Section ### 2 ### :	Cast input fundamental species to mpf arbitary precision datatypes.
	Section ### 3 ### :	Cast input KDs to mpf arbitary precision datatypes.
	Section ### 4 ### :	Define objective function for the minimiser.
	Section ### 5 ### :	Calculate species concentrations from fundamental, and other dependant species.
	Section ### 6 ### :	Return tuple containing the deviation from total fundamental species concentrations.
						More difficult than it first appears as we must count the number of fundamental
						species monomers in each species to know what to multiply the concentration by.
	Section ### 7 ### :	Run the mpmath find_root function on the newly defined objective function, minimising
						the values in returned tuples. This reflects the fundamental species free concentration
						at equilibrium.
	Section ### 8 ### :	Return the dictionary of results, containing concentrations for all species at
						equilibrium.
	
	"""

	assert version_info >= (3, 7), "Requires Python version >=3.7 as dictionary insertion order need to be preserved"
	# Readout denotes the species abundance read out when constructing a BindingSystem
	# object for use in PyBindingCurve simulations/fitting.
	readout = None
	
	# When generated, the custom function string, the custom function, and its
	# arguments are stored as member variables.
	binding_function_string = None
	binding_function=None
	binding_function_arguments=None

	def __init__(self, system_string:str, output_filename=None, dps:int=100):
		"""Construct a minimiser-based custom binding system object

		Args:
			system_string (str): Custom system definition string
			output_filename ([str, Path], optional): Optional file to write generated function to. Defaults to None.
			dps (int, optional): Decimal precision used for calculations performed by MPMath. Defaults to 100.
		"""		
		# Get the reaction dictionary, which takes the form:
		# reaction_dictionary[product]=[[reactant1, reactant2]]
		# list may contain multiple approaches to make product.
		reaction_dict = self.parse_system_definition_string(system_string)
		
		# Now find species (all species present), and fundamental_species
		# dictionaries. They are dictionaries as python >= 3.7 guarantees
		# dicts are ordered. They are used here essentailly like ordered
		# sets, with keys as values for species and values = None.
		species, fundamental_species = self.get_species_and_fundamental_species(reaction_dict)

		# Species_composed_of_matrix is an np array of
		# shape=(len(fundamental_species), len(species)).  Rows and columns
		# ordered by appearance in species and fundamental_species dicts.
		# See function docstring for more complete definition and example.
		species_composed_of_matrix = self.build_species_composed_of_matrix(species, fundamental_species, reaction_dict)
		
		# Get kds needed, in order of appearance in reaction_dict. Again, dict
		# is used to gain functionality of an ordered set.
		kds = self.get_kds(reaction_dict)

		# In cases where there are two routes to make a product, with r1 and
		# r2, or by r3 and r4, the amount must be calculated in a special way:
		# product = r1*r2+r3*r4/(kd1+kd2).
		# Simplified reaction_dict associates one equation with one reaction
		# product for writing out.
		simplified_reaction_dict = self.get_simplified_reaction_dict(reaction_dict, species, fundamental_species)

		# Generate the function string and store in self.custom_func_string
		self.binding_function_string = self.gen_custom_func(species,fundamental_species,kds,simplified_reaction_dict,species_composed_of_matrix)

		# If requested, write the function to a file
		if output_filename is not None:
			out_file = open(output_filename, "w")
			out_file.write("from mpmath import mpf, findroot, mp, almosteq\nmp.dps=100\n\n")
			out_file.write(self.binding_function_string)
			out_file.close()
		
		exec(self.binding_function_string, globals())
		self.binding_function = eval("custom_minimizer_system")
		self.binding_function_arguments = list(signature(self.binding_function).parameters.keys())


	def get_simplified_reaction_dict(self, reaction_dict, species, fundamental_species):
		"""Get simplified reaction dict, required when a product can be made in more than one way

		Args:
			reaction_dict (dict): Reaction dictionary derived from custom binding system string
			species (dict): Dict used as ordered set, containing all species
			fundamental_species (dict): Dict used as ordered set, containing only fundamental species

		Returns:
			dict: Simplified reaction dict, keys are products, values are equations to calculate amounts
		"""		
		srd = {}
		for product, reactions in reaction_dict.items():
			# Only 1 way to make the product (simple case)
			if len(reactions) == 1:
				r1 = reactions[0][0]
				r2 = reactions[0][1]
				if r1 in fundamental_species.keys():
					r1 = f"{r1}_f"
				if r2 in fundamental_species.keys():
					r2 = f"{r2}_f"
				srd[
					product
				] = f"{r1}*{r2}/{self.kd_from_reaction_tuple(reactions[0],product)}"
			else: 
				# 2 or more ways to make the product
				top = ""
				bottom = ""
				for reaction in reactions:
					r1 = reaction[0]
					r2 = reaction[1]
					if r1 in fundamental_species.keys():
						r1 = f"{r1}_f"
					if r2 in fundamental_species.keys():
						r2 = f"{r2}_f"
					top += f"{r1}*{r2}+"
					bottom += f"{self.kd_from_reaction_tuple(reaction,product)}+"
				top = top[:-1]
				bottom = bottom[:-1]
				srd[product] = f"({top})/({bottom})"
		return srd

	def gen_custom_func(self, species: dict, fundamental_species: dict, kds: dict, simplified_reaction_dict: dict, species_composed_of_matrix: np.array):
		"""Generate custom function text

		Args:
			species (dict): All species in system
			fundamental_species (dict): Fundamental species in system
			kds (dict): Unique KDs in order of appearance in the system.
			simplified_reaction_dict (dict): Simplified/unified reaction dictionary
			species_composed_of_matrix (np.array): Numpy int array of monomer counts for all species with shape (len(fundamental_species), len(species))

		Returns:
			str: String representing the custom binding system, solved using MPMath findroot

		The docstring for the MinimizerBindingSystemFactory class outlines the target function text which this function generates.
		Areas denoted by ### 1 ### documented in that text are marked in comments bellow, showing how function construction proceeds.
		"""

		# Section ### 1 ### : Define the custom minimizer function which takes
		# arguments for fundamental species centration and KDs.
		text = (
			f"def custom_minimizer_system("
			+ ",".join(f for f in fundamental_species)
			+ ","
			+ ",".join(kd for kd in kds)
			+ "):\n"
		)
		
		# Section ### 2 ### : Cast input fundamental species to mpf arbitary
		# precision datatypes.
		for fs in fundamental_species:
			text += f"\t{fs}=mpf({fs})\n"
		
		# Section ### 3 ### : Cast input KDs to mpf arbitary precision
		# datatypes.
		for kd in kds:
			text += f"\t{kd}=mpf({kd})\n"
			text +=f"\tif almosteq({kd}, mpf(0), mpf(\"1e-10\")):\n"
			text +=f"\t\t{kd}+=mpf(\"1e-10\")\n"

		
		# Section ### 4 ### : Define objective function for the minimiser.
		text += "\tdef f(" + ",".join(f"{fs}_f" for fs in fundamental_species) + "):\n"
		
		# Section ### 5 ### : Calculate species concentrations from
		# fundamental, and other dependant species.
		for product, reaction in simplified_reaction_dict.items():
			text += f"\t\t{product}={reaction}\n"

		# Section ### 6 ### : Return tuple containing the deviation from total
		# fundamental species concentrations. More difficult than it first
		# appears as we must count the number of fundamental species monomers
		# in each species to know what to multiply the concentration by.
		fundamentals_zero_sum_strs = []
		for fsi, fs in enumerate(fundamental_species):
			monomer_counts = species_composed_of_matrix[:, fsi]
			balance_str = f"{fs}-"
			to_consider = []
			for mci, count in enumerate(monomer_counts):
				if count == 0:
					continue
				monomer_name=list(species.keys())[mci]
				if monomer_name in fundamental_species:
					monomer_name=monomer_name+"_f"
				if count == 1:
					to_consider.append(f"{monomer_name}")
				else:
					to_consider.append(f"{count}*{monomer_name}")
			balance_str += f"({'+'.join(to_consider)})"
			fundamentals_zero_sum_strs.append(balance_str)
		text += "\t\treturn " + ",".join(fundamentals_zero_sum_strs) + "\n"

		# Section ### 7 ### : Run the mpmath find_root function on the newly
		# defined objective function, minimising the values in returned
		# tuples. This reflects the fundamental species free concentration
		# at equilibrium.
		text += (
			"\t"
			+ ",".join(f"{fs}_f" for fs in fundamental_species)
			+ "=findroot(f, ["
			+ ", ".join(f"mpf(0)" for fs in fundamental_species)
			+ "], tol=1e-10, maxsteps=1e6)\n"
		)
		
		# Section ### 5 ### REPEATED: Calculate species concentrations from
		# fundamental, and other dependant species, same as previously written
		# out for objective function
		for product, reaction in simplified_reaction_dict.items():
			text += f"\t{product}={reaction}\n"
		
		# Section ### 8 ### : Return the dictionary of results, containing
		# concentrations for all species at equilibrium.
		text += "\treturn {"
		for fs in fundamental_species:
			text += f"'{fs}_f':{fs}_f,"
		for product in simplified_reaction_dict:
			text += f"'{product}':{product},"
		text = text[:-1] + "}\n"
		
		return text

	def get_kds(self, reaction_dictionary: dict):
		"""Get KDs dictionary from reaction_dictionary

		Args:
				reaction_dictionary (dict): Reaction dictionary

		Returns:
				dict: Dict of KDs written properly. Dict is used as it mimics an ordered set
		"""
		kds = {}
		for product, reactions in reaction_dictionary.items():
			for reaction in reactions:
				kds[self.kd_from_reaction_tuple(reaction, product)] = None
		return kds

	def kd_from_reaction_tuple(self, reaction_tuple: tuple, product: str):
		"""Convenience function to turn (r1, r2), product into kd_r1_r2_product

		Args:
			reaction_tuple (tuple): Reaction tuple, like ('r1', 'r2')
			product (str): product name/symbol

		Returns:
			str: String of the form 'kd_r1_r2_product'
		"""		
		return f"kd_{reaction_tuple[0]}_{reaction_tuple[1]}_{product}"

	def parse_system_definition_string(self, system_string: str):
		"""Parse system definition, generate reaction_dict

		Args:
				system_string (str): Custom system definition
		Returns:
				dict: reaction dictionary product:(reactant1, reactant2)
		"""
		reaction_dict = {}

		lines = [
			l for l in [c.strip() for nl in system_string.lower().split("\n") for c in nl.split(",") if len(c.strip())>=7]
		]
		assert len(lines) > 0, "No system defined"

		# Loop over lines parsing reactions, setting readout, and generating reaction_dict
		for line in lines:
			plus_loc = line.find("+")
			bracket_loc = line.find(">")
			reactant1 = line[0:plus_loc]
			reactant2 = line[plus_loc + 1 : bracket_loc - 2]
			product = line[bracket_loc + 1 :]

			# Check for stars to indicate desired readout
			if "*" in reactant1:
				reactant1 = reactant1.replace("*", "")
				self.readout = reactant1
			if "*" in reactant2:
				reactant2 = reactant2.replace("*", "")
				self.readout = reactant2
			if "*" in product:
				product = product.replace("*", "")
				self.readout = product

			# Add to the objects reaction_dict
			if product not in reaction_dict.keys():
				reaction_dict[product] = []
			reaction_dict[product].append([reactant1, reactant2])

			# Set the readout.
			if self.readout is None:
				first_product = list(reaction_dict.keys())[0]
				print(
					f"* character not found to set readout in custom system, using the first product, which is: {first_product}"
				)
				self.readout = first_product

		return reaction_dict

	def build_species_composed_of_matrix(self, species: dict, fundamental_species: dict, reaction_dict: dict):
		"""Make matrix of fund_species x species containing monomer occurence counts for each species

		Args:
				species (dict): All species dict
				fundamental_species (dict): Fundamental species dict
				reaction_dict (dict): Reaction dictionary

		Returns:
				np.array: Monomer counts array, shape=(len(fundamental_species), len(species))

		# Species_composed_of_matrix is an np array of
		# shape=(len(fundamental_species), len(species)).  Rows and columns
		# ordered by appearance in species and fundamental_species dicts.
		# The following system:
		# p+p<->pp
		# pp+l<->ppl
		# would produce a species_composed_of_matrix as follows:
		#     P L
		#     - - 
		# P  |1 0
		# L  |0 1
		# PP |2 0
		# PPL|2 1
		"""
		len_species = len(species)
		len_fundamental_species = len(fundamental_species)
		species_composed_of_matrix = np.zeros(
			(len_species, len_fundamental_species), dtype=int
		)
		for fsi, _ in enumerate(fundamental_species):
			species_composed_of_matrix[fsi, fsi] = 1

		for si, s in enumerate(list(species.keys())[len_fundamental_species:]):
			species_composed_of_matrix[len_fundamental_species + si, :] = (
				species_composed_of_matrix[
					list(species.keys()).index(reaction_dict[s][0][0]), :
				]
				+ species_composed_of_matrix[
					list(species.keys()).index(reaction_dict[s][0][1]), :
				]
			)
		return species_composed_of_matrix

	def get_species_and_fundamental_species(self, reaction_dict: dict):
		"""Get species and fundamental species dictionaries

		Args:
				reaction_dict (dict): Reaction dictionary

		Returns:
				tuple(dict, dict): Species and fundamental_species dictonaries
		"""
		species = {}
		fundamental_species = {}
		# Populate species
		for product, reactions in reaction_dict.items():
			for reactant1, reactant2 in reactions:
				species[reactant1] = None
				species[reactant2] = None
		for product in reaction_dict:
			species[product] = None

		# Populate fundamental_species
		for s in species.keys():
			if s not in reaction_dict.keys():
				fundamental_species[s] = None
		
		unordered_species=species.copy()
		species={k:v for k, v in fundamental_species.items()}
		species.update({k:v for k,v in unordered_species.items() if k not in fundamental_species.keys()})
		return species, fundamental_species


if __name__ == "__main__":
	custom_system = """
		P+P<->PP
		P+L<->PL
		PP+L<->PPL1
		PP+L<->PPL2
		PPL1+L<->PPL1L2
		PPL2+L<->PPL1L2
		"""
	custom_generated_binding_system = MinimizerBindingSystemFactory(
		custom_system, output_filename="tmp.py"
	)
	print(custom_generated_binding_system.binding_function_arguments)
	print(custom_generated_binding_system.binding_function(10,10,10,10,10,10,10,10))
