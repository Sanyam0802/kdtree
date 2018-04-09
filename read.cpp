#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

typedef vector<vector<double>> points;

points readData()
{
	ifstream dataFile("dataset.txt");
	string line;
	int dimension = 0;
	int numberOfPoints = 0;
	points all_points;


	if (dataFile.is_open())
    {
    	int i = 0;

    	while (! dataFile.eof() )
    	{
    		stringstream ss;
      		getline (dataFile,line);
      		ss.str(line);

      		if(line == ""){
      			continue;
      		}

    		if (i==0){
    			ss>>dimension;
    			ss>>numberOfPoints;
    			i++;
    			continue;
    		}     		
    		vector<double> point;
    		for(int j=0;j<dimension;j++){
    			double x;
    			ss>>x;
    			point.push_back(x);
    		}
    		all_points.push_back(point);
     		i++;
    	}
    }

    return all_points;
}

int main(){
	points p = readData();
	cout<<p.size()<<endl;
	return 0;
}