#include "Delaunay.h"
#include <ctime>
int main()
{
	clock_t start, end;
	double total_time;

	std::string file_name = "leaf";
	int n = 10000;
	Delaunay M(file_name, n);

	start = clock();

	M.Triangulation();

	end = clock();
	total_time = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << " run time: " << total_time << std::endl;
	// M.~Delaunay();
	return 0;
}
