#include "Analyzer.h"
#include <csignal>

bool do_break;
void KeyboardInterrupt_endJob(int signum) {
    do_break = true;
}

void usage() {
  cout << endl;
  cout << "You have entered too little arguments, please type:\n";
  cout << "./Analyzer infile.root outfile.root\n";
  cout << "or\n";
  cout << "./Analyzer -out outfile.root -inlist infile.root infile.root infile.root\n";
  cout << "or\n";
  cout << "./Analyzer -out outfile.root -in infile.root\n";
  cout << "Available options are:\n";
  cout << "-CR: to run over the control regions (not the usual output)\n";
  cout << "-C: use a different config folder than \n";
  cout << "\n";
  cout << "Note the -inlist option must be used with the -out option and will add all files, which are after the flag (breaks on '-').";

  exit(EXIT_FAILURE);
}

void parseCommandLine(int argc, char *argv[], vector<string> &inputnames, string &outputname, bool &setCR, string &configFolder) {
  if (argc==1)
    usage();
  if(argc < 3)
    usage();
  for (int arg=1; arg<argc; arg++) {

    if (strcmp(argv[arg], "-CR") == 0) {
      setCR = true;
      arg++;
      continue;
    }else if (strcmp(argv[arg], "-C") == 0) {
      configFolder=argv[arg+1];
      cout << "Analyser: ConfigFolder " << configFolder << endl;
      arg++;
      continue;
    }else if (strcmp(argv[arg], "-in") == 0) {
      inputnames.push_back(argv[arg+1]);
      cout << "Analyser: Inputfile " << inputnames.back() << endl;
      arg++;
      continue;
    }else if (strcmp(argv[arg], "-out") == 0) {
      outputname=argv[arg+1];
      cout << "Analyser: Outputfile " << outputname << endl;
      arg++;
      continue;
    }else if (strcmp(argv[arg], "-inlist") == 0) {
      inputnames.push_back(argv[arg+1]);
      cout << "Analyser: Inputfilelist " << inputnames.back() << endl;
      arg++;
      while( arg<argc and (argv[arg][0] != '-')){
        inputnames.push_back(argv[arg]);
        cout << "Analyser: Inputfilelist " << inputnames.back() << endl;
        arg++;
      }
      if(argv[arg][0] == '-')
        arg--;
      continue;
    }else if(inputnames.size()==0){
      inputnames.push_back(argv[arg]);
    }else if(outputname==""){
      outputname = argv[arg];
    }
  }

  ifstream ifile(inputnames.back());
  if ( !ifile && inputnames.back().find("root://") == string::npos && inputnames.back().find("root\\://") == string::npos) {
    std::cout << "The file '" << inputnames.back() << "' doesn't exist" << std::endl;
    exit(EXIT_FAILURE);
  }
  return;
}

int main (int argc, char* argv[]) {

  bool setCR = false;
  do_break =false;

  string outputname;
  string configFolder="PartDet";
  vector<string> inputnames;


  //get the command line options in a nice loop
  parseCommandLine(argc, argv, inputnames, outputname, setCR, configFolder);

  //setup the analyser
  Analyzer testing(inputnames, outputname, setCR, configFolder);

  //catch ctrl+c and just exit the loop
  //this way we still have the output
  signal(SIGINT,KeyboardInterrupt_endJob);


  for(int i=0; i < testing.nentries; i++) {
    testing.clear_values();
    testing.preprocess(i);
    testing.fill_histogram();
    //this will be set if ctrl+c is pressed
    if(do_break){
      testing.nentries=i;
      break;
    }
  }
  testing.printCuts();
  return 0;
}
