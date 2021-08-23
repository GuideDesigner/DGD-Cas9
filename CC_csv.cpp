#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <set>

class ccClass{
private:
public:
  ccClass(){
    id="";
    pos_a=0;
    pos_b=0;
    cluster_id=0;
  };
  ~ccClass(){};

  std::string id;
  int pos_a;
  int pos_b;
  int cluster_id;
  int cluster_number;
  
  void print(){
    std::cout << id << ","
	      << pos_a << ","
	      << pos_b << ","
	      << cluster_id << ","
	      << cluster_number << "\n";
	      
  }
};



std::set < std::pair < int , int > > find_seed ( std::map < int , std::map < int , ccClass > > container ){
  int pos_a;
  int pos_b;
  std::pair < int , int > seed;
  std::set < std::pair < int , int > > seed_set;
  
  for ( auto i : container ){
    pos_a = i.first;
    for ( auto j : i.second ){
      pos_b = j.first;
      if ( container[pos_a-1][pos_b+1].id.size() == 0 ){
	seed.first = pos_a;
	seed.second = pos_b;
	seed_set.insert(seed);
      }
    }
  }
  return seed_set;
};


void explore ( std::map < int , std::map < int , ccClass > > & container,
	       std::set < std::pair < int , int > > seed_set){
  int pos_a;
  int pos_b;
  
  int id_number;
  for ( auto i : seed_set ){
    pos_a = i.first;
    pos_b = i.second;
    id_number = container[pos_a][pos_b].cluster_id;
    int a = pos_a + 1;
    int b = pos_b - 1;

    while ( container[a][b].id.size() > 0 ){
      container[a][b].cluster_id = id_number;
      a ++;
      b --;
    }
    

  }
}

void scan_CC(std::map < int , std::map < int , ccClass > > & container, int & id_number){
  std::set < std::pair < int , int > > seed_set = find_seed ( container );
  int pos_a;
  int pos_b;

  // designate id
  for ( auto i : seed_set ){
    pos_a = i.first;
    pos_b = i.second;
    container[pos_a][pos_b].cluster_id = id_number;
    id_number ++;
  }

  explore(container, seed_set);
};

int main ( int argc , char ** argv )
{
  if ( argc == 1 ){
    std::cout << "Usage CC_csv input.csv\n";
    return 0;
  }
  // std::string::size_type sz;
  std::string infile(argv[1]);
  
  std::string id;
  int pos_a;
  int pos_b;
  
  std::string tmp;
  
  std::map < int , std::map < int , ccClass > > fav;
  std::map < int , std::map < int , ccClass > > dis;

  std::ifstream fin(infile.c_str());

  std::getline ( fin , tmp );

  while ( std::getline ( fin , tmp , ',' ) ){
    id = tmp;
    std::getline ( fin , tmp , ',' );
    pos_a = atoi ( tmp.c_str() );
    std::getline ( fin , tmp );
    pos_b = atoi ( tmp.c_str() );

    fav[pos_a][pos_b].id=id;
    fav[pos_a][pos_b].pos_a=pos_a;
    fav[pos_a][pos_b].pos_b=pos_b;
    
  }
         
  fin.close();

  int cluster_id = 1;
  scan_CC(fav,cluster_id);
  scan_CC(dis,cluster_id);

  std::map < int , std::map < int , ccClass > > annotated;
		   
  std::map < int , int > CC_number_map;

  for ( auto i : fav ){
    for ( auto j : i.second ){
      if ( j.second.cluster_id > 0 ) {
	annotated[i.first][j.first]=j.second;
	CC_number_map[j.second.cluster_id]++;
      }
    }
  }
  for ( auto i : dis ){
    for ( auto j : i.second ){
      if ( j.second.cluster_id > 0 ) {
	annotated[i.first][j.first]=j.second;
	CC_number_map[j.second.cluster_id]++;
      }
    }
  }

  for ( auto i : annotated ){
    for ( auto j : i.second ){
      int ccid=annotated[i.first][j.first].cluster_id;
      annotated[i.first][j.first].cluster_number = CC_number_map[ccid];
    }
  }
  
  fin.open(infile.c_str());
  std::getline ( fin , tmp );
  std::cout << tmp << ",CC_id,CC_num\n" ;

  while ( std::getline ( fin , tmp , ',' ) ){
    std::getline ( fin , tmp , ',' );
    pos_a = atoi ( tmp.c_str() );
    std::getline ( fin , tmp );
    pos_b = atoi ( tmp.c_str() );

    annotated[pos_a][pos_b].id=id;
    annotated[pos_a][pos_b].pos_a=pos_a;
    annotated[pos_a][pos_b].pos_b=pos_b;
    annotated[pos_a][pos_b].print();
  }

  fin.close();

  return 0;
}
