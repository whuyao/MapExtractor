#pragma warning(disable:4068)

#include "QtCore"


#include <omp.h>
#include <iostream>
#include <mapnik/map.hpp>
#include <mapnik/load_map.hpp>
#include <mapnik/layer.hpp>
#include <mapnik/rule.hpp>
#include <mapnik/feature_type_style.hpp>
#include <mapnik/symbolizer.hpp>
#include <mapnik/text/placements/dummy.hpp>
#include <mapnik/text/text_properties.hpp>
#include <mapnik/text/formatting/text.hpp>
#include <mapnik/datasource_cache.hpp>
#include <mapnik/font_engine_freetype.hpp>
#include <mapnik/agg_renderer.hpp>
#include <mapnik/expression.hpp>
#include <mapnik/color_factory.hpp>
#include <mapnik/image_util.hpp>
#include <mapnik/unicode.hpp>
#include <mapnik/save_map.hpp>
#include <mapnik/datasource_cache.hpp>
#include <mapnik/cairo_io.hpp>
#include <mapnik/projection.hpp>
#include <mapnik/geometry.hpp>
#include <windows.h>
#include <string>
#include "GeoHash.h"

#if defined(HAVE_CAIRO)
#include <mapnik/cairo/cairo_renderer.hpp>
#include <mapnik/cairo/cairo_image_util.hpp>
#endif

using namespace mapnik;


int oneLayerPerLoc()
{
	try
	{
		// =============== set parameters from console ===============
		//if (argc < 11)
		//{
		// 		std::cout << "please input parameters like this: " << std::endl;
		// 		std::cout << ">> OutputTile.exe MapnikXMLFile LocationFile SaveImageDir ImageRadius PerLocImgNum StartLocIndex RandomDis HashPrecision ThreadsNum" << std::endl;
		// 		std::cout << "exit." << std::endl;
		// 		return -1;
		//}
		//const char* sMapXmlFile = argv[1];
		//const char* sLatLonFile = argv[2];
		//const char* sSaveImageDir = argv[3];
		//int nImgSize = atoi(argv[4]);
		//double dSquareDis = atof(argv[5]);
		//int nPerLocImgNum = atoi(argv[6]);
		//int nStartFID = atoi(argv[7]);
		//double dRandomDis = atoi(argv[8]);
		//int nHashPrescision = atoi(argv[9]);
		//int nThreadsNum = atoi(argv[10]);
		// ==============================================

		// =============== set parameters ===============
		const char* sMapXmlFile = "osm_xml/mapnik.xml";
		const char* sLatLonFile = "./data/lonlat.csv";
		const char* sSaveImageDir = "./osm_image_test";
		int nImgSize = 256;
		double dSquareDis = 100;
		int nPerLocImgNum = 20;
		int nStartFID = 0;
		double dRandomDis = dSquareDis * 0.1;
		int nHashPrescision = 12;
		int nThreadsNum = 6;
		// ==============================================

		datasource_cache::instance().register_datasources("./mapnik/input");

		const std::string srs_merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs +over";
		//const std::string srs_merc_rotate = "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0";
		std::string srs_rotate_proj_base = "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=";
		const std::string srs_lonlat = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

		QList<Map> mvMaps;

		std::cout << "loading map from xml ..." << std::endl;
		Map _map(256, 256);
		load_map(_map, sMapXmlFile);
		std::cout << "load map finish." << std::endl;

		std::cout << "Querying from postgis ..." << std::endl << std::endl;

		for (int i = 0; i < nThreadsNum; i++)
			mvMaps.append(_map);

		QFile _infile(sLatLonFile);
		if (!_infile.open(QIODevice::ReadOnly))
		{
			std::cout << "open latlon file error!" << std::endl;
			return false;
		}
		QTextStream _in(&_infile);
		QList<QString> svAllPoints;

		_in.readLine();

		while (!_in.atEnd())
		{
			QString smsg = _in.readLine();
			svAllPoints.append(smsg);
		}

#pragma omp parallel for num_threads(nThreadsNum) schedule(dynamic)
		for (int k = nStartFID; k < svAllPoints.size(); k++)
		{
			int nThreadNum = omp_get_thread_num();
			box2d<double> square_bounds(0, 0, dSquareDis * 2, dSquareDis * 2);

			QString smsg = svAllPoints[k];
			QStringList smsglist = smsg.split(",");

			int nFID = smsglist[0].trimmed().toInt();
			double dx = smsglist[1].trimmed().toDouble();
			double dy = smsglist[2].trimmed().toDouble();

			std::string strLatLngHash = GeoHash::EncodeLatLon(dy, dx, nHashPrescision);

			geometry::point<double> dpoint(dx, dy);

			QString sDirPath = QString("%1/%2").arg(sSaveImageDir).arg(strLatLngHash.data());
			QDir _dir;
			if (!_dir.exists(sDirPath))
				_dir.mkdir(sDirPath);

			srand((unsigned)time(NULL));

			for (int i = 0; i < nPerLocImgNum; i++)
			{

				if (i == 0)
				{
					mvMaps[nThreadNum].set_srs(srs_merc);
					projection merc = projection(srs_merc);
					projection lonlat = projection(srs_lonlat);
					proj_transform tr(lonlat, merc);
					geometry::point<double> dTmpPoint(dpoint);
					tr.forward(dTmpPoint);
					square_bounds.re_center(dTmpPoint.x, dTmpPoint.y);
				}
				else
				{
					double dRandLon = (rand() / double(RAND_MAX)) * 360;
					std::string srs_rotate_proj = srs_rotate_proj_base;
					std::stringstream _strstream;
					_strstream << dRandLon;
					srs_rotate_proj += _strstream.str();

					double dRandX = (rand() / double(RAND_MAX)) * dRandomDis;
					double dRandY = (rand() / double(RAND_MAX)) * dRandomDis;


					mvMaps[nThreadNum].set_srs(srs_rotate_proj);
					geometry::point<double> dTmpPoint(dpoint);
					projection rotate_merc = projection(srs_rotate_proj);
					projection lonlat = projection(srs_lonlat);
					proj_transform tr(lonlat, rotate_merc);
					tr.forward(dTmpPoint);
					square_bounds.re_center(dTmpPoint.x + dRandX, dTmpPoint.y + dRandY);

				}

				mvMaps[nThreadNum].zoom_to_box(square_bounds);

				mapnik::image_rgba8 buf(mvMaps[nThreadNum].width(), mvMaps[nThreadNum].height());
				mapnik::agg_renderer<mapnik::image_rgba8> r(mvMaps[nThreadNum], buf);
				r.apply();

				QString sImgPath = QString("%1/%2_%3.jpg").arg(sDirPath).arg(strLatLngHash.data()).arg(i);
				const char* sImgFileName = sImgPath.toStdString().data();

				std::cout << "saving image of FID " << nFID << " : " << sImgFileName << " ..." << std::endl;

				try
				{
					mapnik::save_to_file(buf, sImgFileName, "jpeg");
					std::cout << "save success." << std::endl << std::endl;
				}
				catch (std::exception const& ex)
				{
					std::cerr << "### std::exception: " << ex.what();
					std::cout << " " << sImgFileName << "\nnow retrying ..." << std::endl;
					i--;
					continue;
				}

			}
		}
	}
	catch (std::exception const& ex)
	{
		std::cerr << "### std::exception: " << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...)
	{
		std::cerr << "### Unknown exception." << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
	//将多层xml加载然后进行渲染

	try
	{
		// =============== set parameters from console ===============
// 		if (argc < 11)
// 		{
// 			std::cout << "please input parameters like this: " << std::endl;
// 			std::cout << ">> OutputTile.exe MapnikXMLDir LocationFile SaveImageDir ImageSize ImageRadius PerLocImgNum StartLocIndex RandomDis HashPrecision ThreadsNum" << std::endl;
// 			std::cout << "exit." << std::endl;
// 			return -1;
// 		}
// 		const char* sMapXmlDir = argv[1];
// 		const char* sLatLonFile = argv[2];
// 		const char* sSaveImageDir = argv[3];
// 		int nImgSize = atoi(argv[4]);
// 		double dSquareDis = atof(argv[5]);
// 		int nPerLocImgNum = atoi(argv[6]);
// 		int nStartFID = atoi(argv[7]);
// 		double dRandomDis = atoi(argv[8]);
// 		int nHashPrescision = atoi(argv[9]);
// 		int nThreadsNum = atoi(argv[10]);
		// ==============================================

		// =============== set parameters ===============
		const char* sMapXmlDir = "./osm_xml";            //xml文件夹
		const char* sLatLonFile = "./data/lonlat.csv";   //经纬度列表文件
		const char* sSaveImageDir = "./osm_image_test";  //保存图片结果文件夹
		int nImgSize = 256;                              //输出图片大小
		double dSquareDis = 100;                         //输出范围半径大小（m）
		int nPerLocImgNum = 20;                          //每个位置随机输出多少张图片
		int nStartFID = 0;                               //程序从哪个FID的经纬度开始运行
		double dRandomDis = dSquareDis * 0.1;            //随机采样时移动的范围大小
		int nHashPrescision = 12;                        //输出文件夹和文件名的GeoHash精度
		int nThreadsNum = 6;                             //使用多少个线程
		// ==============================================

		//注册数据源
		datasource_cache::instance().register_datasources("./mapnik/input");

		//墨卡托投影
		const std::string srs_merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs +over";
		
		//用于旋转的投影的基本部分，后面加上0~360即可进行随机旋转
		std::string srs_rotate_proj_base = "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=";

		//WGS84经纬度坐标
		const std::string srs_lonlat = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

		QList<Map> mvMaps;
		QList<QList<Map>> mapvAllMaps;

		//遍历xml文件加载所有的xml
		QDir _dir(sMapXmlDir);
		QFileInfoList finfo = _dir.entryInfoList(QStringList(QString("*.xml")), QDir::Files, QDir::Name);
		if (finfo.isEmpty())
		{
			std::cout << "not find file error" << std::endl;
			return EXIT_FAILURE;
		}

		int nXMLFileNum = finfo.size();
		std::cout << "XML file num = " << finfo.size() << std::endl;
		std::cout << "loading map from xml ..." << std::endl;
		foreach(QFileInfo f1, finfo)
		{
			std::cout << "loading " << f1.baseName().toStdString() << std::endl;
			QString sXmlFilePath = f1.absoluteFilePath();
			std::string strMapXmlFile = sXmlFilePath.toStdString();
			Map _map(nImgSize, nImgSize);
			load_map(_map, strMapXmlFile);
			mvMaps.append(_map);
		}
		std::cout << "load map finish." << std::endl;
		std::cout << "Querying from postgis ..." << std::endl << std::endl;

		//为了使用多线程，新建nThreadsNum个Map列表
		for (int i = 0; i < nThreadsNum; i++)
			mapvAllMaps.append(mvMaps);

		//读取经纬度列表
		QFile _infile(sLatLonFile);
		if (!_infile.open(QIODevice::ReadOnly))
		{
			std::cout << "open location file error!" << std::endl;
			return EXIT_FAILURE;
		}
		QTextStream _in(&_infile);
		QList<QString> svAllPoints;

		_in.readLine();

		while (!_in.atEnd())
		{
			QString smsg = _in.readLine();
			svAllPoints.append(smsg);
		}

#pragma omp parallel for num_threads(nThreadsNum) schedule(dynamic)
		for (int k = nStartFID; k < svAllPoints.size(); k++)
		{
			int nThreadID = omp_get_thread_num();
			box2d<double> square_bounds(0, 0, dSquareDis * 2, dSquareDis * 2);

			QString smsg = svAllPoints[k];
			QStringList smsglist = smsg.split(",");

			int nFID = smsglist[0].trimmed().toInt();
			double dx = smsglist[1].trimmed().toDouble();
			double dy = smsglist[2].trimmed().toDouble();

			std::string strLatLngHash = GeoHash::EncodeLatLon(dy, dx, nHashPrescision);

			geometry::point<double> dpoint(dx, dy);

			QString sDirPath = QString("%1/%2").arg(sSaveImageDir).arg(strLatLngHash.data());
			QDir _dir;
			if (!_dir.exists(sDirPath))
				_dir.mkdir(sDirPath);

			srand((unsigned)time(NULL));

			for (int i = 0; i < nPerLocImgNum; i++)
			{
				std::string srs_proj;
				double dRandX = 0, dRandY = 0;

				if (i == 0)
					srs_proj = srs_merc;
				else
				{
					double dRandLon = (rand() / double(RAND_MAX)) * 360;
					srs_proj = srs_rotate_proj_base;
					std::stringstream _strstream;
					_strstream << dRandLon;
					srs_proj += _strstream.str();

					dRandX = (rand() / double(RAND_MAX)) * dRandomDis;
					dRandY = (rand() / double(RAND_MAX)) * dRandomDis;
				}

				projection merc = projection(srs_proj);
				projection lonlat = projection(srs_lonlat);
				proj_transform tr(lonlat, merc);
				geometry::point<double> dTmpPoint(dpoint);
				tr.forward(dTmpPoint);
				square_bounds.re_center(dTmpPoint.x + dRandX, dTmpPoint.y + dRandY);

				for (int j = 0; j < nXMLFileNum; j++)
				{
					mapvAllMaps[nThreadID][j].set_srs(srs_proj);
					mapvAllMaps[nThreadID][j].zoom_to_box(square_bounds);

					mapnik::image_rgba8 buf(mapvAllMaps[nThreadID][j].width(), mapvAllMaps[nThreadID][j].height());
					mapnik::agg_renderer<mapnik::image_rgba8> r(mapvAllMaps[nThreadID][j], buf);
					r.apply();

					QString sImgPath = QString("%1/%2_%3_%4.jpg").arg(sDirPath).arg(strLatLngHash.data()).arg(i).arg(j);
					const char* sImgFileName = sImgPath.toStdString().data();

					std::cout << "saving image of FID " << nFID << " : " << sImgFileName << " ..." << std::endl;

					try
					{
						mapnik::save_to_file(buf, sImgFileName, "jpeg");
						std::cout << "save success." << std::endl << std::endl;
					}
					catch (std::exception const& ex)
					{
						std::cerr << "### std::exception: " << ex.what();
						std::cout << " " << sImgFileName << "\nnow retrying ..." << std::endl;
						j--;
						continue;
					}

				}
			}
		}
	}
	catch (std::exception const& ex)
	{
		std::cerr << "### std::exception: " << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...)
	{
		std::cerr << "### Unknown exception." << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
