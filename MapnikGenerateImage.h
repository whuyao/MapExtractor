#ifndef MAPNIK_GENERATE_IMAGE_H_FILE
#define MAPNIK_GENERATE_IMAGE_H_FILE

#pragma warning(disable:4068)

#include <mapnik/map.hpp>
#include <QtCore>

bool GenerateImage(mapnik::Map mapOut, const char* sOutFileName);

bool LoadXML2Map(const char* sMapXmlDir, QList<mapnik::Map>& mvMaps, int nImgSize);

bool LoadLocation(const char* sLocationFile, QList<QString>& svAllPoints);

bool GeneratePerLocOneImage(
	const char* sDataSourceDir,//数据源存放位置，为mapnik库的lib文件夹中的input文件夹
	const char* sMapXmlDir,    //xml文件夹
	const char* sLocationFile, //经纬度列表文件
	const char* sSaveImageDir, //保存图片结果文件夹
	int nImageSize,            //输出的图片大小
	double dSquareDis,         //输出范围半径大小(m)
	int nStartFID,             //程序从哪个FID的经纬度开始运行
	int nHashPrecision,        //输出文件名的GeoHash精度
	int nThreadsNum            //使用多少个线程进行输出
);

bool GeneratePerLocMultiImage(
	const char* sDataSourceDir,
	const char* sMapXmlDir,
	const char* sLocationFile,
	const char* sSaveImageDir,
	int nImageSize,
	double dSquareDis,
	int nStartFID,
	int nHashPrecision,
	int nThreadsNum,
	int nRandImageNum
);

// 当图片过多时，每次训练载入时都需要很长时间，所以将所有图片的文件名列在文件中
bool OutputObjectList(const char* sImageDir, const char* sObjectListFile);

// 检查输出的图片是否有缺漏
bool CheckPerLocOneImageResult(const char* sSaveImageDir, const char* sLocationFile, int nImageLayer, int nHashPrecision, const char* sMissingLocFile);

// 检查输出的图片是否有缺漏
bool CheckPerLocMultiImageResult(const char* sSaveImageDir, const char* sLocationFile, int nImageLayer, int nHashPrecision, int nRandImageNum, const char* sMissingLocFile);

bool CheckPerLocMultiImageResult(const char* sLabelFile, const char* sObjectListFile, const char* sMissingFile);


// 按顺序生成label
bool outputLabelFile(const char* sLocFile, const char* sLabelFile);


#endif
