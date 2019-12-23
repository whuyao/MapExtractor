#pragma warning(disable:4068)

#include <QtCore>
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

#include "MapnikGenerateImage.h"

#if defined(HAVE_CAIRO)
#include <mapnik/cairo/cairo_renderer.hpp>
#include <mapnik/cairo/cairo_image_util.hpp>
#endif

using namespace mapnik;


int main(int argc, char *argv[])
{
	const char* sDataSourceDir = "./mapnik/input";
	const char* sMapXmlDir = "./origin_osm_xml"; // "./osm_xml";            //xml文件夹
	const char* sLocationFile = "./data/wuhan.csv";   //经纬度列表文件
	const char* sSaveImageDir = "./data/osm_per_loc_multi_image_wuhan"; // "./data/osm_per_loc_multi_image_wuhan";  //保存图片结果文件夹
	int nImgSize = 256;                              //输出图片大小
	double dSquareDis = 100;                         //输出范围半径大小（m）
	int nStartFID = 0;                               //程序从哪个FID的经纬度开始运行
	int nHashPrecision = 12;                        //输出文件夹和文件名的GeoHash精度
	int nThreadsNum = 6;                             //使用多少个线程
	int nRandImageNum = 20;                          //每个位置随机生成多少张图
	int nImageLayer = 8;                             //输出图片的层数，即xml文件的数量

	// if directory is not exist then create it
	QDir dir(sSaveImageDir);
	if (!dir.exists())
		dir.mkpath(QFileInfo(sSaveImageDir).absoluteFilePath());

	//return 1;

	// ================= generate images =====================
	//GeneratePerLocOneImage(sDataSourceDir, sMapXmlDir, sLocationFile, sSaveImageDir, nImgSize, dSquareDis, nStartFID, nHashPrecision, nThreadsNum);
	GeneratePerLocMultiImage(sDataSourceDir, sMapXmlDir, sLocationFile, sSaveImageDir, nImgSize, dSquareDis, nStartFID, nHashPrecision, nThreadsNum, nRandImageNum);


	// ================= check the output result =====================
	const char* sMissingLocFile = "./data/osm_iamges_missing.csv";
	//CheckPerLocOneImageResult(sSaveImageDir, sLocationFile, nImageLayer, nHashPrecision, sMissingLocFile);
	CheckPerLocMultiImageResult(sSaveImageDir, sLocationFile, nImageLayer, nHashPrecision, nRandImageNum, sMissingLocFile);


	// ================= output the object list of iamges =====================
	const char* sObjectListFile = "./data/osm_image_test_object_list.txt";
	OutputObjectList(sSaveImageDir, sObjectListFile);


	// ================= output the labels of iamges =====================
	const char* sLabelFile = "./data/osm_iamge_label.csv";
	outputLabelFile(sLocationFile, sLabelFile);

	

	return 0;
}
