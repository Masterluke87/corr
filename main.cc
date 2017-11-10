#include <fstream>
#include <iostream>
#include <sstream>
#include "ops_io.h"

void read_input(std::ifstream* inputfile, std::string* sysfile, int nroe,int llim, int ulim, std::string* wavefile)
{
  std::string tmpline;
  std::stringstream ss_input;
  if (inputfile->is_open())
  {
    while(std::getline((*inputfile),tmpline))
    {
      //remove comments
      if (tmpline.find_first_of("#")!=std::string::npos)
        tmpline.erase(tmpline.find_first_of("#"),tmpline.length());
      if (tmpline.find_first_of(" ")!=std::string::npos)
        tmpline.erase(tmpline.find_first_of(" "),tmpline.length());
      if (tmpline.compare("\n")!=0)
        ss_input <<tmpline<<" ";
    }
  }
  inputfile->close();
  ss_input >> (*sysfile) >> nroe >> llim >> ulim >> (*wavefile);
  std::cout<<"Sysfile    : "<<(*sysfile)<<std::endl;
  std::cout<<"#Electrons : "<<nroe<<std::endl;
  std::cout<<"llim       : "<<llim<<std::endl;
  std::cout<<"ulim       : "<<ulim<<std::endl;
  std::cout<<"Wavefile   : "<<(*wavefile)<<std::endl;
  std::cout<<"=======================================\n\n";



}


int main(int argc, char const *argv[]) {
  int     nroe;           //Nr of electrons  (if negative, read in center of mass)
  int     llim;           //first MO used for correlation
  int     ulim;           //last  MO used for correlation
  int    nroao;
  int     nroa;
  long long  int nrofint;      //Nr of two electron Integrals

  std::string sysfile;
  std::string wavefile;
  std::cout<<"PSI4/MP2 Program \n";

  if(argc != 3){
    std::cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  std::ifstream inputfile (argv[1]);

  read_input(&inputfile,&sysfile,nroe,llim,ulim,&wavefile);
  get_sys_size(sysfile, &nroao, &nroa,  &nrofint);

  std::cout<< "System sizes read from       : " << sysfile << "\n";
  std::cout<< "Nr of basis functions        : " << nroao   << "\n";
  std::cout<< "Nr of atoms                  : " << nroa << "\n";
  std::cout<< "Nr of non zero 2el integrals : " << nrofint << "\n";
  std::cout<< "=======================================\n\n";

  std::cout << "--END--" << '\n';
  return 0;
}
