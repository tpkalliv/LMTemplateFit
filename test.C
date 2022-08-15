void outputs (TString *arrNames);


void test() {

	TString arrNames[] = {"yes", "no"};


	outputs(arrNames);

}

void outputs (TString *arrNames) {

	cout << "1" << endl;
	cout << arrNames[0] << endl;

	/*
	for (int i = 0; i < sizeof(values); i++) values[i] = 1 + i;

	
	for (int j = 0; j < sizeof(names); j++) cout << "Name " << j << " is " << names[j] << endl;
	for (int k = 0; k < sizeof(values); k++) cout << "Value " << k << " is " << values[k] << endl;	
	*/
}