# Project: connection_table

CPP  = g++
CC   = gcc
OBJ  = atom.o emptyBin.o find_no_fflines.o frame_read_xmolout.o get_settings.o identify_atom_types.o identify_bond_parameters.o string_conversion_functions.o find_no_types_input_atoms.o make_connection_table.o find_valence_params.o identify_molecules.o identify_molecule_name.o identify_ch3.o identify_c3.o identify_double_bonds.o identify_header_lines.o connection_table.o read_table.o read_reax_table.o identify_carbon_3_neighbors.o find_no_connt_lines.o identify_atoms_in_molecule.o write_reaction_input_file.o bin_atoms.o ring_analysis.o identify_voids.o aux_functions.o $(RES)
LINKOBJ  = atom.o emptyBin.o find_no_fflines.o frame_read_xmolout.o get_settings.o identify_atom_types.o identify_bond_parameters.o string_conversion_functions.o find_no_types_input_atoms.o make_connection_table.o find_valence_params.o identify_molecules.o identify_molecule_name.o identify_ch3.o identify_c3.o identify_double_bonds.o identify_header_lines.o connection_table.o read_table.o read_reax_table.o identify_carbon_3_neighbors.o find_no_connt_lines.o identify_atoms_in_molecule.o write_reaction_input_file.o bin_atoms.o ring_analysis.o identify_voids.o aux_functions.o $(RES)
LIBS =  -L /usr/lib  
INCS =  -I /usr/include/ 
#CXXINCS = -I /usr/include/c++/3.4.6/backward -I /usr/include/c++/3.4.2 -g 
CXXINCS = -I /usr/include/c++/3.4.6/backward -I /usr/include/c++/3.4.2 -std=c++11
BIN  = connection_table
CXXFLAGS = $(CXXINCS)  
CFLAGS = $(INCS)  
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before connection_table all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o "debug_connection_table" $(LIBS)

atom.o: atom.cpp
	$(CPP) -c atom.cpp -o atom.o $(CXXFLAGS)

emptyBin.o: emptyBin.cpp
	$(CPP) -c emptyBin.cpp -o emptyBin.o $(CXXFLAGS)

find_no_fflines.o: find_no_fflines.cpp
	$(CPP) -c find_no_fflines.cpp -o find_no_fflines.o $(CXXFLAGS)

frame_read_xmolout.o: frame_read_xmolout.cpp
	$(CPP) -c frame_read_xmolout.cpp -o frame_read_xmolout.o $(CXXFLAGS)

get_settings.o: get_settings.cpp
	$(CPP) -c get_settings.cpp -o get_settings.o $(CXXFLAGS)

identify_atom_types.o: identify_atom_types.cpp
	$(CPP) -c identify_atom_types.cpp -o identify_atom_types.o $(CXXFLAGS)

identify_bond_parameters.o: identify_bond_parameters.cpp
	$(CPP) -c identify_bond_parameters.cpp -o identify_bond_parameters.o $(CXXFLAGS)

string_conversion_functions.o: string_conversion_functions.cpp
	$(CPP) -c string_conversion_functions.cpp -o string_conversion_functions.o $(CXXFLAGS)

find_no_types_input_atoms.o: find_no_types_input_atoms.cpp
	$(CPP) -c find_no_types_input_atoms.cpp -o find_no_types_input_atoms.o $(CXXFLAGS)

make_connection_table.o: make_connection_table.cpp
	$(CPP) -c make_connection_table.cpp -o make_connection_table.o $(CXXFLAGS)

find_valence_params.o: find_valence_params.cpp
	$(CPP) -c find_valence_params.cpp -o find_valence_params.o $(CXXFLAGS)

identify_molecules.o: identify_molecules.cpp
	$(CPP) -c identify_molecules.cpp -o identify_molecules.o $(CXXFLAGS)

identify_molecule_name.o: identify_molecule_name.cpp
	$(CPP) -c identify_molecule_name.cpp -o identify_molecule_name.o $(CXXFLAGS)

identify_ch3.o: identify_ch3.cpp
	$(CPP) -c identify_ch3.cpp -o identify_ch3.o $(CXXFLAGS)

identify_c3.o: identify_c3.cpp
	$(CPP) -c identify_c3.cpp -o identify_c3.o $(CXXFLAGS)

identify_double_bonds.o: identify_double_bonds.cpp
	$(CPP) -c identify_double_bonds.cpp -o identify_double_bonds.o $(CXXFLAGS)

identify_header_lines.o: identify_header_lines.cpp
	$(CPP) -c identify_header_lines.cpp -o identify_header_lines.o $(CXXFLAGS)

connection_table.o: connection_table.cpp
	$(CPP) -c connection_table.cpp -o connection_table.o $(CXXFLAGS)

read_table.o: read_table.cpp
	$(CPP) -c read_table.cpp -o read_table.o $(CXXFLAGS)

read_reax_table.o: read_reax_table.cpp
	$(CPP) -c read_reax_table.cpp -o read_reax_table.o $(CXXFLAGS)

identify_carbon_3_neighbors.o: identify_carbon_3_neighbors.cpp
	$(CPP) -c identify_carbon_3_neighbors.cpp -o identify_carbon_3_neighbors.o $(CXXFLAGS)

find_no_connt_lines.o: find_no_connt_lines.cpp
	$(CPP) -c find_no_connt_lines.cpp -o find_no_connt_lines.o $(CXXFLAGS)

identify_atoms_in_molecule.o: identify_atoms_in_molecule.cpp
	$(CPP) -c identify_atoms_in_molecule.cpp -o identify_atoms_in_molecule.o $(CXXFLAGS)

write_reaction_input_file.o: write_reaction_input_file.cpp
	$(CPP) -c write_reaction_input_file.cpp -o write_reaction_input_file.o $(CXXFLAGS)

bin_atoms.o: bin_atoms.cpp
	$(CPP) -c bin_atoms.cpp -o bin_atoms.o $(CXXFLAGS)

ring_analysis.o: ring_analysis.cpp
	$(CPP) -c ring_analysis.cpp -o ring_analysis.o $(CXXFLAGS)

identify_voids.o: identify_voids.cpp
	$(CPP) -c identify_voids.cpp -o identify_voids.o $(CXXFLAGS)

aux_functions.o: aux_functions.cpp
	$(CPP) -c aux_functions.cpp -o aux_functions.o $(CXXFLAGS)
