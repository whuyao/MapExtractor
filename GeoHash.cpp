#include "GeoHash.h"
#include <algorithm>

#define PI       3.14159265358979323846   // pi

const string Base32 = "0123456789bcdefghjkmnpqrstuvwxyz";
const int Bits[5] = { 16, 8, 4, 2, 1 };
const string Neighbors[2][4] = {
	{
		"p0r21436x8zb9dcf5h7kjnmqesgutwvy", // Top
		"bc01fg45238967deuvhjyznpkmstqrwx", // Right
		"14365h7k9dcfesgujnmqp0r2twvyx8zb", // Bottom
		"238967debc01fg45kmstqrwxuvhjyznp", // Left
	},
	{
		"bc01fg45238967deuvhjyznpkmstqrwx", // Top
		"p0r21436x8zb9dcf5h7kjnmqesgutwvy", // Right
		"238967debc01fg45kmstqrwxuvhjyznp", // Bottom
		"14365h7k9dcfesgujnmqp0r2twvyx8zb", // Left
	}
};

const string Borders[2][4] = {
	{ "prxz", "bcfguvyz", "028b", "0145hjnp" },
	{ "bcfguvyz", "prxz", "0145hjnp", "028b" }
};

#define EARTH_RADIUS 6378.137


//角度转弧度
double rad(double d)
{
	return d * PI / 180.0;
}

/*******************************************************************
* 返回两个地点P1(lng1, lat1)和P2(lng2, lat2)的近似大地线距离(单位 km)
* double lng1: 第一个地点的经度(角度)
* double lat1: 第一个地点的纬度(角度)
* double lng2: 第二个地点的经度(角度)
* double lat2: 第二个地点的纬度(角度)
* 返回值: double, 两个点之间的近似大地线距离，单位 km
*******************************************************************/
double getPtsDist(double lat1, double lng1, double lat2, double lng2)
{
	double radLat1 = rad(lat1);
	double radLat2 = rad(lat2);
	double radLon1 = rad(lng1);
	double radLon2 = rad(lng2);
	double a = radLat1 - radLat2;
	double b = radLon1 - radLon2;

	double s = 2 * asinf(sqrtf(sin(a / 2)*sin(a / 2) + cos(radLat1)*cos(radLat2)*sin(b / 2)*sin(b / 2)));
	s = s * EARTH_RADIUS;

	return s;
}

std::string GeoHash::CalculateAdjacent(string hash, Direction direction)
{
	//hash = hash;
	transform(hash.begin(), hash.end(), hash.begin(), ::tolower);

	char lastChr = hash[hash.length() - 1];
	int type = hash.length() % 2;
	int dir = (int)direction;
	string nHash = hash.substr(0, hash.length() - 1);

	if (Borders[type][dir].find(lastChr) != -1)
	{
		nHash = CalculateAdjacent(nHash, (Direction)dir);
	}
	return nHash + Base32[Neighbors[type][dir].find(lastChr)];
}

void GeoHash::RefineInterval(double* interval, int cd, int mask)
{
	if ((cd & mask) != 0)
	{
		interval[0] = (interval[0] + interval[1]) / 2;
	}
	else
	{
		interval[1] = (interval[0] + interval[1]) / 2;
	}
}

void GeoHash::DecodeGeoHash(string geohash, double& dlat, double& dlon)
{
	bool even = true;
	double lat[2] = { -90.0, 90.0 };
	double lon[2] = { -180.0, 180.0 };

	for (char c : geohash)
	{
		int cd = Base32.find(c);
		for (int j = 0; j < 5; j++)
		{
			int mask = Bits[j];
			if (even)
			{
				RefineInterval(lon, cd, mask);
			}
			else
			{
				RefineInterval(lat, cd, mask);
			}
			even = !even;
		}
	}

	dlat = (lat[0] + lat[1]) / 2;
	dlon = (lon[0] + lon[1]) / 2;
}

std::string GeoHash::EncodeLatLon(double latitude, double longitude, int precision /*= 12*/)
{
	bool even = true;
	int bit = 0;
	int ch = 0;
	string geohash = "";

	double lat[2] = { -90.0, 90.0 };
	double lon[2] = { -180.0, 180.0 };

	if (precision < 1 || precision > 20) precision = 12;

	while (geohash.length() < precision)
	{
		double mid;

		if (even)
		{
			mid = (lon[0] + lon[1]) / 2;
			if (longitude > mid)
			{
				ch |= Bits[bit];
				lon[0] = mid;
			}
			else
				lon[1] = mid;
		}
		else
		{
			mid = (lat[0] + lat[1]) / 2;
			if (latitude > mid)
			{
				ch |= Bits[bit];
				lat[0] = mid;
			}
			else
				lat[1] = mid;
		}

		even = !even;
		if (bit < 4)
			bit++;
		else
		{
			geohash += Base32[ch];
			bit = 0;
			ch = 0;
		}
	}
	return geohash;
}
	
double GeoHash::distanceBetweenGeoHashes(string hash1, string hash2)
{
	double lat1, lon1, lat2, lon2;
	DecodeGeoHash(hash1, lat1, lon1);
	DecodeGeoHash(hash2, lat2, lon2);
	return getPtsDist(lat1, lon1, lat2, lon2);
}

