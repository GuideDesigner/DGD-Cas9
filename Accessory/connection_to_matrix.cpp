#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

int main ( int argc , char ** argv ){
  if ( argc == 1 ) {
    std::cout << "Usage:  connection_to_matrix connection_file.txt sequence_size >  outfile\n";
    return 0;
  }

  std::string infile=argv[1];
  int seq_size = atoi(argv[2]);
  std::string tmp;
  int value;
  int x,y;

  std::ifstream fin ( infile );
  
  std::vector < std::string > header;
  // std::map < int , std::map < int , double > > val_mat;

  /**
   * Reading header
   */
  
  fin >> tmp ;
  for ( int i = 1 ; i < seq_size ; i ++ ){
    fin >> tmp ;
    header.push_back(tmp);
  }

  /**
   * Writing header
   */
  std::cout << "ID";
  for ( int i = 1 ; i < seq_size ; i ++ ){
    for ( int j = 1 ; j < seq_size ; j ++ ){
      std::cout << "," << "Connection_" << header[i] << "_" << header[j];
    }
  }
  std::cout << "\n";


  /**
   * Reading and writing matrix
   */

  while ( fin >> tmp ){
    std::map < int , std::map < int , int > > val_mat;
    for ( int i = 1 ; i < seq_size ; i ++ ){
      fin >> value;
      x = i + 1 ;
      y = value ;
      val_mat[x][y] += 1;
    }
    std::cout << tmp ;
    for ( int i = 1 ; i < seq_size ; i ++ ){
      for ( int j = 1 ; j < seq_size ; j ++ ){
     	std::cout << "," << val_mat[i+1][j+1];
      }
    }
    std::cout << "\n";
  }
  
  
  return 0;
}
