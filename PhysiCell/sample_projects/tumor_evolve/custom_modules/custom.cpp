/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

Cell_Definition TC;

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 

	//birth behavior
	
	if(parameters.strings("enable_sigmoid") == "1"){
	  cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_sigmoid;
	}else if(parameters.strings("enable_carryingcapacity") == "1"){
	  cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_carryingcapacity;
	 }else{
	  cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_standard;
	}


	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 

	*/ 

	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0;

	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.functions.cycle_model = live;
	if(parameters.doubles("cell_migration_speed") > 0.0){
	  cell_defaults.phenotype.motility.is_motile = true;
	}else{
	  cell_defaults.phenotype.motility.is_motile = false;
	}

	cell_defaults.custom_data.add_variable( "growth_mut" , "dimensionless", 0 );
	//Division value
	cell_defaults.custom_data.add_variable( "growth" , "dimensionless", 0 );
	//Deleterious mutation barcode 
	cell_defaults.custom_data.add_variable( "mig_mut" , "dimensionless", 0 );
	//migration rate
	cell_defaults.custom_data.add_variable( "mig" , "dimensionless", parameters.doubles("cell_migration_speed" ));
	//Stores the cell's current pressure
	cell_defaults.custom_data.add_variable( "gf" , "dimensionless", 0 );

	TC = cell_defaults;
	TC.custom_data[0] = 0;
	TC.custom_data[1] = TC.phenotype.cycle.data.transition_rate( 0, 0 ); 
	TC.custom_data[2] = 0;
	TC.custom_data[3] = TC.phenotype.motility.migration_speed;
	TC.custom_data[4] = 0;
	
        cell_defaults.type = 0;
	TC.type = 1;
	TC.parameters.pReference_live_phenotype = &( TC.phenotype );

	TC.phenotype.cycle.data.transition_rate( 0, 0 ) *= 1; //return here
	TC.custom_data[1] = TC.phenotype.cycle.data.transition_rate( 0, 0 ); 
	TC.phenotype.death.rates[apoptosis_model_index] = parameters.doubles("cell_death_rate") * TC.phenotype.cycle.data.transition_rate( 0, 0 ); 

	TC.phenotype.motility.migration_speed = parameters.doubles("cell_migration_speed" );
	TC.custom_data[3] = TC.phenotype.motility.migration_speed;
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 

  //default_microenvironment_options.X_range = {-500, 500}; 
  //default_microenvironment_options.Y_range = {-500, 500};
  //default_microenvironment_options.Z_range = {0, 0};
	  
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	  Cell* pC;

	  double x = 0.0;                                                      
  //  for ( x = 0; x < 1; x = x + 1 ) {                                 
	  for ( x = 0; x < parameters.doubles("initial_cell_number"); x = x + 1 ) {

	    pC = create_cell(TC );                    
	    pC->assign_position( 0.0 + NormalRandom(0, 20), 0.0 + NormalRandom(0, 20), 0.0 );
	    
	  }

	
	return; 
}

void update_cell_and_death_parameters_standard( Cell* pCell, Phenotype& phenotype, double dt )
{

  if( phenotype.death.dead == true )
    { return; }

  double pres = pCell->state.simple_pressure;
  pCell->custom_data[4] = pres;

  if( pres > parameters.doubles( "pressure_threshold" ) )
    {
      pCell->phenotype.cycle.data.transition_rate( 0, 0 ) = 0.0;
    }else{
      pCell->phenotype.cycle.data.transition_rate( 0, 0 ) = pCell->custom_data[1];
  }

  return;
}

//This function is to smooth the cut off between 0 and full birth rates as a function of 
//pressure. This function encodes a sigmoidal relationship between pressure and birth
void update_cell_and_death_parameters_sigmoid( Cell* pCell, Phenotype& phenotype, double dt )
{

  if( phenotype.death.dead == true )
    { return; }

  double pres = pCell->state.simple_pressure;
  pCell->custom_data[4] = pres;
  pCell->phenotype.cycle.data.transition_rate( 0, 0 ) = pCell->custom_data[1] * (1 - 1/(1 +  exp(-5 * (pres - 1))));

  return;
  
}

//This function updates the death rate as the tumor grows so its growth ultimately slows and
//reaches a carrying capacity
void update_cell_and_death_parameters_carryingcapacity( Cell* pCell, Phenotype& phenotype, double dt )
{

  if( phenotype.death.dead == true )
    { return; }

  int apoptosis_model_index = pCell->phenotype.death.find_death_model_index( "Apoptosis" );

  double f_carrying_capacity =  (*all_cells).size() / parameters.doubles("cellNumEndCondition") ;

  //  std::cout<< "carrying capacity: " + std::to_string(f_carrying_capacity) + "\n";
  
  double relative_death = parameters.doubles("cell_death_rate") +  0.75 * (1 - parameters.doubles("cell_death_rate"))/( 1 + exp(-10 * (f_carrying_capacity - 0.5)));

    pCell->phenotype.death.rates[apoptosis_model_index] = relative_death * pCell->custom_data[1];
    
    double pres = pCell->state.simple_pressure;
    pCell->custom_data[4] = pres;

    if( pres > parameters.doubles( "pressure_threshold" ) )
      {
	pCell->phenotype.cycle.data.transition_rate( 0, 0 ) = 0.0;
      }else{
        pCell->phenotype.cycle.data.transition_rate( 0, 0 ) = pCell->custom_data[1];
    }

  return;
}



bool stopCondition()
{

  return((*all_cells).size() <= parameters.doubles("cellNumEndCondition") & (*all_cells).size() > 0);

}


std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
