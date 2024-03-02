//
// Created by vivek on 14/11/18.
//

#include "sharpness.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    // get the command line arguments
	if (argc != 3){
        cout<<"Three arguments expected"<<endl;
        return 1;
    }

	float scale_factor = atof(argv[3]);

	assess_sharpness(argv[1], argv[2], scale_factor);

    return 0;
}
