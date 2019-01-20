#ifndef settings_H
#define settings_H
#include <string>

using std::string;
class settings{
      public:
             int table_flag;
             int molecule_flag;
             int c3_flag;
             int ring_flag;
             int void_flag;
             int saturate_H_flag;
             int periodic_flag;
             double grid_size;
             double bo_cut_off;

         settings(){
             molecule_flag = 0; //0 means no analysis
             table_flag = c3_flag = ring_flag = void_flag = saturate_H_flag = periodic_flag = 0;
             grid_size = 4.0;
             bo_cut_off = 0.3;
         }
             
         int return_molecule_flag(){ return molecule_flag;}
         int return_c3_flag(){ return c3_flag;}
         int return_ring_flag(){ return ring_flag;}
         int return_void_flag(){ return void_flag;}
         int return_saturate_H_flag(){ return saturate_H_flag;}
         int return_periodic_flag(){ return periodic_flag;}
         double return_grid_size(){ return grid_size;}
         
};
      
#endif
