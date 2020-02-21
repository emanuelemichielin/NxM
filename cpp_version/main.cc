#include "header.h"
#include "TApplication.h"





void StandaloneApplication(int argc, char** argv) {
    Driver();
}


int main(int argc, char** argv) {

    if(argc >1){
	if(std::atoi(argv[1]) == 0){
	    Driver();
	    return 0;
	}
    }else{
        TApplication app("ROOT Application", &argc, argv);
        StandaloneApplication(app.Argc(), app.Argv());
        app.Run();
        Driver();
        return 0;
    }
}

