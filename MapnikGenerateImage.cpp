
#include "MapnikGenerateImage.h"
#include <mapnik/agg_renderer.hpp>
#include <mapnik/save_map.hpp>
#include <mapnik/image_util.hpp>
#include <mapnik/load_map.hpp>
#include <mapnik/projection.hpp>
#include <mapnik/datasource_cache.hpp>
#include "GeoHash.h"
#include <omp.h>
using namespace mapnik;

bool GenerateImage(mapnik::Map mapOut, const char* sOutFileName)
{
	mapnik::image_rgba8 buf(mapOut.width(), mapOut.height());
	mapnik::agg_renderer<mapnik::image_rgba8> r(mapOut, buf);
	r.apply();

	std::cout << "saving image of " << sOutFileName << " ... " << std::endl;

	try
	{
		mapnik::save_to_file(buf, sOutFileName, "jpeg");
		std::cout << "save success." << std::endl << std::endl;
	}
	catch (std::exception const& ex)
	{
		std::cerr << "### std::exception: " << ex.what();
		std::cout << "\nnow retrying ... " << std::endl;
		return false;
	}
	
	return true;
}

bool LoadXML2Map(const char* sMapXmlDir, QList<mapnik::Map>& mvMaps, int nImgSize)
{
	QDir _dir(sMapXmlDir);
	QFileInfoList finfo = _dir.entryInfoList(QStringList(QString("*.xml")), QDir::Files, QDir::Name);
	if (finfo.isEmpty())
	{
		std::cout << "Not find xml file in " << sMapXmlDir << std::endl;
		return false;
	}

	int nXmlFileNum = finfo.size();
	std::cout << "XML file num = " << finfo.size() << std::endl;
	std::cout << "Loading map from xml ..." << std::endl;

	foreach(QFileInfo f1, finfo)
	{
		std::cout << "loading " << f1.baseName().toStdString() << std::endl;
		QString sXmlFilePath = f1.absoluteFilePath();
		std::string strMapXmlFile = sXmlFilePath.toStdString();
		mapnik::Map _map(nImgSize, nImgSize);
		mapnik::load_map(_map, strMapXmlFile);
		mvMaps.append(_map);
	}
	std::cout << "load map finish." << std::endl;

	return true;
}

bool LoadLocation(const char* sLocationFile, QList<QString>& svAllPoints)
{
	QFile _infile(sLocationFile);
	if (!_infile.open(QIODevice::ReadOnly))
	{
		std::cout << "open location file error!" << std::endl;
		return false;
	}
	QTextStream _in(&_infile);
	
	_in.readLine();
	while (!_in.atEnd())
	{
		QString smsg = _in.readLine();
		svAllPoints.append(smsg);
	}

	return true;
}

bool GeneratePerLocOneImage(const char* sDataSourceDir, const char* sMapXmlDir, const char* sLocationFile, const char* sSaveImageDir,  int nImageSize,  double dSquareDis, int nStartFID,  int nHashPrecision,int nThreadsNum)
{
	try
	{
		//注册数据源
		datasource_cache::instance().register_datasources(sDataSourceDir);
		//墨卡托投影
		const std::string srs_merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs +over";
		//WGS84经纬度坐标
		const std::string srs_lonlat = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

		QList<mapnik::Map> mvMaps;
		QList<QList<mapnik::Map>> mapvAllMaps;

		if (!LoadXML2Map(sMapXmlDir, mvMaps, nImageSize))
		{
			return false;
		}

		QList<QString> svAllPoints;
		if (!LoadLocation(sLocationFile, svAllPoints))
		{
			return false;
		}

		QDir _outdir;
		if (!_outdir.exists(sSaveImageDir))
		{
			std::cout << "The iamge save dir " << sSaveImageDir << " did not exists, now create it." << std::endl;
			_outdir.mkdir(sSaveImageDir);
		}

		//为了使用多线程，新建nThreadsNum个Map列表
		for (int i = 0; i < nThreadsNum; i++)
			mapvAllMaps.append(mvMaps);

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

			std::string strLatLngHash = GeoHash::EncodeLatLon(dy, dx, nHashPrecision);

			geometry::point<double> dpoint(dx, dy);

			projection merc = projection(srs_merc);
			projection lonlat = projection(srs_lonlat);
			proj_transform tr(lonlat, merc);
			tr.forward(dpoint);
			square_bounds.re_center(dpoint.x, dpoint.y);

			for (int j = 0; j < mapvAllMaps[0].size(); j++)
			{
				mapvAllMaps[nThreadID][j].set_srs(srs_merc);
				mapvAllMaps[nThreadID][j].zoom_to_box(square_bounds);

				mapnik::image_rgba8 buf(mapvAllMaps[nThreadID][j].width(), mapvAllMaps[nThreadID][j].height());
				mapnik::agg_renderer<mapnik::image_rgba8> r(mapvAllMaps[nThreadID][j], buf);
				r.apply();

				QString sImgPath = QString("%1/%2_%3.jpg").arg(sSaveImageDir).arg(strLatLngHash.data()).arg(j);
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
	catch (std::exception const& ex)
	{
		std::cerr << "### std::exception: " << ex.what() << std::endl;
		return false;
	}
	catch (...)
	{
		std::cerr << "### Unknown exception." << std::endl;
		return false;
	}

	return true;
}

bool GeneratePerLocMultiImage(const char* sDataSourceDir, const char* sMapXmlDir, const char* sLocationFile, const char* sSaveImageDir, int nImageSize, double dSquareDis, int nStartFID, int nHashPrecision, int nThreadsNum, int nRandImageNum)
{
	try
	{
		//注册数据源
		datasource_cache::instance().register_datasources(sDataSourceDir);
		//墨卡托投影
		std::string srs_merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs +over";
		//WGS84经纬度坐标
		const std::string srs_lonlat = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
		//旋转坐标系的前半段
		std::string srs_rotate_proj_base = "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=";
		double dRandomDis = dSquareDis * 0.1;

		QList<mapnik::Map> mvMaps;
		QList<QList<mapnik::Map>> mapvAllMaps;

		if (!LoadXML2Map(sMapXmlDir, mvMaps, nImageSize))
		{
			return false;
		}

		QList<QString> svAllPoints;
		if (!LoadLocation(sLocationFile, svAllPoints))
		{
			return false;
		}

		QDir _outdir;
		if (!_outdir.exists(sSaveImageDir))
		{
			std::cout << "The iamge save dir " << sSaveImageDir << " did not exists, now create it." << std::endl;
			_outdir.mkdir(sSaveImageDir);
		}

		//为了使用多线程，新建nThreadsNum个Map列表
		for (int i = 0; i < nThreadsNum; i++)
			mapvAllMaps.append(mvMaps);

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

			std::string strLatLngHash = GeoHash::EncodeLatLon(dy, dx, nHashPrecision);

			for (int i = 0; i < nRandImageNum; i++)
			{
				std::string srs_proj = srs_merc;
				geometry::point<double> dpoint(dx, dy);

				double dRandX = 0;
				double dRandY = 0;

				if (i != 0)
				{
					double dRandLon = (rand() / double(RAND_MAX)) * 360;
					dRandX = (rand() / double(RAND_MAX)) * dRandomDis;
					dRandY = (rand() / double(RAND_MAX)) * dRandomDis;

					srs_proj = srs_rotate_proj_base;
					std::stringstream _strstream;
					_strstream << dRandLon;
					srs_proj += _strstream.str();
				}

				projection merc = projection(srs_proj);
				projection lonlat = projection(srs_lonlat);
				proj_transform tr(lonlat, merc);
				tr.forward(dpoint);
				square_bounds.re_center(dpoint.x + dRandX, dpoint.y + dRandY);

				for (int j = 0; j < mapvAllMaps[0].size(); j++)
				{
					mapvAllMaps[nThreadID][j].set_srs(srs_proj);
					mapvAllMaps[nThreadID][j].zoom_to_box(square_bounds);

					mapnik::image_rgba8 buf(mapvAllMaps[nThreadID][j].width(), mapvAllMaps[nThreadID][j].height());
					mapnik::agg_renderer<mapnik::image_rgba8> r(mapvAllMaps[nThreadID][j], buf);
					r.apply();

					QString sImgPath = QString("%1/%2_%3_%4.jpg").arg(sSaveImageDir).arg(strLatLngHash.data()).arg(i).arg(j);
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
		return false;
	}
	catch (...)
	{
		std::cerr << "### Unknown exception." << std::endl;
		return false;
	}

	return true;
}

bool OutputObjectList(const char* sImageDir, const char* sObjectListFile)
{
	QDir _dir(sImageDir);
	QFileInfoList finfolist = _dir.entryInfoList(QStringList(QString("*.jpg")), QDir::Files, QDir::Name);

	if (finfolist.isEmpty())
	{
		std::cout << "not find files in " << sImageDir << std::endl;
		return false;
	}

	QFile _outfile(sObjectListFile);
	if (!_outfile.open(QIODevice::WriteOnly))
	{
		std::cout << "Create object list file error!" << std::endl;
		return false;
	}
	QTextStream _out(&_outfile);

	foreach(QFileInfo finfo, finfolist)
	{
		QString sFileName = finfo.fileName();
		//std::cout << sFileName.toStdString() << std::endl;

		_out << sFileName << "\r\n";
	}
	_out.flush();

	_outfile.close();
	return true;
}

bool CheckPerLocOneImageResult(const char* sSaveImageDir, const char* sLocationFile, int nImageLayer, int nHashPrecision, const char* sMissingLocFile)
{
	QList<QString> svLocGeoHash;
	QFile _infile(sLocationFile);
	if (!_infile.open(QIODevice::ReadOnly))
	{
		std::cout << "open location file error!" << std::endl;
		return false;
	}
	QTextStream _in(&_infile);

	QFile _outfile(sMissingLocFile);
	if (!_outfile.open(QIODevice::WriteOnly))
	{
		std::cout << "Create Missing Location File error!" << std::endl;
		return false;
	}
	QTextStream _out(&_outfile);

	QString sHeadLine = _in.readLine();
	_out << sHeadLine << "\r\n";

	while (!_in.atEnd())
	{
		QString smsg = _in.readLine();
		QStringList smsglist = smsg.split(",");
		double dx = smsglist[1].trimmed().toDouble();
		double dy = smsglist[2].trimmed().toDouble();

		std::string strLatLngHash = GeoHash::EncodeLatLon(dy, dx, nHashPrecision);
		QString _hash = QString("%1").arg(strLatLngHash.data());

		svLocGeoHash.push_back(_hash);
	}

	QHash<QString, std::vector<std::string>> hash_objects;
	QDir _dir(sSaveImageDir);
	QFileInfoList finfolist = _dir.entryInfoList(QStringList(QString("*.jpg")), QDir::Files, QDir::Name);
	if (finfolist.isEmpty())
	{
		std::cout << "Not find files in " << sSaveImageDir << std::endl;
		return false;
	}
	foreach(QFileInfo finfo, finfolist)
	{
		std::string sFilePath = finfo.filePath().toStdString();
		QString sBaseName = finfo.baseName();
		QStringList sBaseNameList = sBaseName.split("_");
		QString sGeoHash = sBaseNameList[0];

		if (hash_objects.contains(sGeoHash))
		{
			hash_objects[sGeoHash].push_back(sFilePath);
		}
		else
		{
			std::vector<std::string> svPerLocFilePath;
			svPerLocFilePath.push_back(sFilePath);
			hash_objects.insert(sGeoHash, svPerLocFilePath);
		}
	}

	int _count = 0;
	for (int i = 0; i < svLocGeoHash.size(); i++)
	{
		QString sGeoHash = svLocGeoHash[i];

		if ((!hash_objects.contains(sGeoHash)) || (hash_objects[sGeoHash].size() != nImageLayer))
		{
			std::cout << "Missing images in " << sGeoHash.toStdString() << std::endl;
			double dlng = 0, dlat = 0;
			GeoHash::DecodeGeoHash(sGeoHash.toStdString(), dlat, dlng);
			QString _outstr = QString("%1,%2,%3").arg(_count).arg(dlng, 0, 'f', 8).arg(dlat, 0, 'f', 8);
			_out << _outstr << "\r\n";
			_count++;
		}
	}

	_outfile.close();
	_infile.close();
	return true;
}

bool CheckPerLocMultiImageResult(const char* sSaveImageDir, const char* sLocationFile, int nImageLayer, int nHashPrecision, int nRandImageNum, const char* sMissingLocFile)
{
	QList<QString> svLocGeoHash;
	QFile _infile(sLocationFile);
	if (!_infile.open(QIODevice::ReadOnly))
	{
		std::cout << "open location file error!" << std::endl;
		return false;
	}
	QTextStream _in(&_infile);

	QFile _outfile(sMissingLocFile);
	if (!_outfile.open(QIODevice::WriteOnly))
	{
		std::cout << "Create Missing Location File error!" << std::endl;
		return false;
	}
	QTextStream _out(&_outfile);

	QString sHeadLine = _in.readLine();
	_out << sHeadLine << "\r\n";

	while (!_in.atEnd())
	{
		QString smsg = _in.readLine();
		QStringList smsglist = smsg.split(",");
		double dx = smsglist[1].trimmed().toDouble();
		double dy = smsglist[2].trimmed().toDouble();

		std::string strLatLngHash = GeoHash::EncodeLatLon(dy, dx, nHashPrecision);
		QString _hash = QString("%1").arg(strLatLngHash.data());

		svLocGeoHash.push_back(_hash);
	}

	QHash<QString, std::vector<std::string>> hash_objects;
	QDir _dir(sSaveImageDir);
	QFileInfoList finfolist = _dir.entryInfoList(QStringList(QString("*.jpg")), QDir::Files, QDir::Name);
	if (finfolist.isEmpty())
	{
		std::cout << "Not find files in " << sSaveImageDir << std::endl;
		return false;
	}
	foreach(QFileInfo finfo, finfolist)
	{
		std::string sFilePath = finfo.filePath().toStdString();
		QString sBaseName = finfo.baseName();
		QStringList sBaseNameList = sBaseName.split("_");
		QString sGeoHash = sBaseNameList[0];

		if (hash_objects.contains(sGeoHash))
		{
			hash_objects[sGeoHash].push_back(sFilePath);
		}
		else
		{
			std::vector<std::string> svPerLocFilePath;
			svPerLocFilePath.push_back(sFilePath);
			hash_objects.insert(sGeoHash, svPerLocFilePath);
		}
	}

	int _count = 0;
	for (int i = 0; i < svLocGeoHash.size(); i++)
	{
		QString sGeoHash = svLocGeoHash[i];

		if ((!hash_objects.contains(sGeoHash)) || (hash_objects[sGeoHash].size() != nImageLayer * nRandImageNum))
		{
			std::cout << "Missing images in " << sGeoHash.toStdString() << std::endl;
			double dlng = 0, dlat = 0;
			GeoHash::DecodeGeoHash(sGeoHash.toStdString(), dlat, dlng);
			QString _outstr = QString("%1,%2,%3").arg(_count).arg(dlng, 0, 'f', 8).arg(dlat, 0, 'f', 8);
			_out << _outstr << "\r\n";
			_count++;
		}
	}

	_outfile.close();
	_infile.close();

	return true;
}

bool CheckPerLocMultiImageResult(const char* sLabelFile, const char* sObjectListFile, const char* sMissingFile)
{
	QFile _labelfile(sLabelFile);
	if (!_labelfile.open(QIODevice::ReadOnly))
	{
		std::cout << "Open label file error!" << std::endl;
		return false;
	}
	QTextStream _inlabel(&_labelfile);

	QFile _ojfile(sObjectListFile);
	if (!_ojfile.open(QIODevice::ReadOnly))
	{
		std::cout << "open object list file error!" << std::endl;
		return false;
	}
	QTextStream _inoj(&_ojfile);

	QList<QString> svLabelHashList;

	_inlabel.readLine();
	while (!_inlabel.atEnd())
	{
		QString smsg = _inlabel.readLine();
		QStringList smsglist = smsg.split(",");

		QString sGeoHash = smsglist[0].trimmed();
		int nLabel = smsglist[1].trimmed().toInt();

		if (svLabelHashList.contains(sGeoHash))
		{
			std::cout << sGeoHash.toStdString() << ", " << nLabel << std::endl;
			std::cout << svLabelHashList.indexOf(sGeoHash) << std::endl;
		}
		else
			svLabelHashList.push_back(sGeoHash);
	}

	_ojfile.close();
	_labelfile.close();
	return true;
}

bool outputLabelFile(const char* sLocFile, const char* sLabelFile)
{
	QFile _infile(sLocFile);
	if (!_infile.open(QIODevice::ReadOnly))
	{
		std::cout << "open location file error!" << std::endl;
		return false;
	}
	QTextStream _in(&_infile);

	QFile _outfile(sLabelFile);
	if (!_outfile.open(QIODevice::WriteOnly))
	{
		std::cout << "create label file error!" << std::endl;
		return false;
	}
	QTextStream _out(&_outfile);

	_in.readLine();
	_out << "GeoHash,Label\r\n";

	int _count = 0;
	while (!_in.atEnd())
	{
		QString smsg = _in.readLine();
		QStringList smsglist = smsg.split(",");

		double dlng = smsglist[1].trimmed().toDouble();
		double dlat = smsglist[2].trimmed().toDouble();

		std::string strLatLngHash = GeoHash::EncodeLatLon(dlat, dlng, 12);

		_out << strLatLngHash.data() << "," << _count << "\r\n";
		_count++;
	}


	_outfile.close();
	_infile.close();
	return true;
}

