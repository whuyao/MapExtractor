#ifndef GEOHASH_WRITE_BY_YAO_VER_10
#define GEOHASH_WRITE_BY_YAO_VER_10

#include <string>
using namespace std;

/*******************************************************************
* 返回两个地点P1(lng1, lat1)和P2(lng2, lat2)的近似大地线距离(单位 km)
* double lng1: 第一个地点的经度(角度)
* double lat1: 第一个地点的纬度(角度)
* double lng2: 第二个地点的经度(角度)
* double lat2: 第二个地点的纬度(角度)
* 返回值: double, 两个点之间的近似大地线距离，单位 km
*******************************************************************/
double getPtsDist(double lat1, double lng1, double lat2, double lng2);

namespace GeoHash {
	enum Direction 
	{
		Top = 0,
		Right = 1,
		Bottom = 2,
		Left = 3,
		Limit = 4
	};
		
	void RefineInterval(double* interval, int cd, int mask);

	// reference: https://www.cnblogs.com/aiweixiao/p/6188081.html
	string CalculateAdjacent(string hash, Direction direction);
	void DecodeGeoHash(string geohash, double& dlat, double& dlon);
	string EncodeLatLon(double latitude, double longitude, int precision = 12);

	// return unit: KM
	double distanceBetweenGeoHashes(string hash1, string hash2);
};




#endif




