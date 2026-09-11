#pragma once

// #include <Ponca/Ponca>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

template <typename Point>
std::vector<Point> ReadFile(const std::string& fileName)
{
    int N          = 0;
    int D          = 0;
    int hasNormal = false;

    std::ifstream file(fileName);

    if (!file.is_open())
        throw std::runtime_error("Can not open: " + fileName);

    std::string header;
    if (!std::getline(file, header))
        throw std::runtime_error("Empty file");

    std::istringstream headerStream(header);
    headerStream >> N >> D >> hasNormal;

    std::cout << header << std::endl;
    std::cout << N << ", " << D << ", " << std::boolalpha << hasNormal << std::endl;
    if (headerStream.bad())
        throw std::runtime_error("Invalid header format");

    if (D != Point::Dim)
        throw std::runtime_error("Wrong point dimension. Found: " + std::to_string(D) +
                                 ", expected: " + std::to_string(Point::Dim));

    if (N <= 0)
        throw std::runtime_error("The number of point must be positive.");

    std::vector<Point> points;
    points.resize(N);

    for (int i = 0; i < N; ++i)
    {
        std::string line;
        if (!std::getline(file, line))
            throw std::runtime_error("Unexpected end of file.");

        auto& pos    = points[i].pos();
        auto& normal = points[i].normal();
        std::istringstream lineStream(line);

        for (int j = 0; j < D; ++j)
            if (!(lineStream >> pos[j]))
                throw std::runtime_error("Missing coordinate or invalid format at line " + std::to_string(i));

        if (hasNormal)
            for (int j = 0; j < D; ++j)
                if (!(lineStream >> normal[j]))
                    throw std::runtime_error("Missing coordinate or invalid format at line " + std::to_string(i));
    }

    return points;
}

/** WRITE JSON FORMAT **/

template <typename Vector>
void WriteVector(std::ofstream& file, const Vector& v)
{
    if (v.rows() == 0)
        return;

    file << v[0];
    for (unsigned int i = 1; i < v.rows(); ++i)
        file << "," << v[i];
}

template <typename Point>
void WriteArray(std::ofstream& file, const char* name, const std::vector<typename Point::Scalar>& values)
{
    if (values.empty())
        return;

    file << '"' << name << "\":{\"N\":" << values.size() << ",\"D\": 1,\"data\":[" << values[0];
    for (unsigned int i = 1; i < values.size(); ++i)
        file << ',' << values[i];
    file << "]}";
}

template <typename Point>
void WriteArray(std::ofstream& file, const char* name, const std::vector<typename Point::VectorType>& values)
{
    if (values.empty())
        return;

    file << '"' << name << "\":{\"N\":" << values.size() << ",\"D\":" << Point::Dim << ",\"data\":["    ;
    WriteVector(file, values[0]);
    for (unsigned int i = 1; i < values.size(); ++i)
    {
        file << ',';
        WriteVector(file, values[i]);
    }
    file << "]}";
}


template <typename Point, bool Normals>
void WriteArray(std::ofstream& file, const char* name, const std::vector<Point>& points)
{
    file << "\"" << name << "\":{\"N\":" << points.size() << ",\"D\":" << Point::Dim << ",";
    {
        file << "\"pos\": [";
        for (size_t i = 0; i < points.size(); ++i)
        {   
            for (size_t j = 0; j < Point::Dim; ++j)
            {
                if (i != 0 || j != 0) file << ",";
                file << points[i].pos()[j];
            }
        }
        file << "]";
    }

    if constexpr (Normals)
    {
        file << ",\"normals\": [";
        for (size_t i = 0; i < points.size(); ++i)
        {   
            for (size_t j = 0; j < Point::Dim; ++j)
            {
                if (i != 0 || j != 0) file << ",";
                file << points[i].normal()[j];
            }
        }
        file << "]";
    }
    file << "}";
}

template<typename Point, typename Input, typename Output>
void WriteResult(std::ofstream& file, const char* method, const char* function, const std::vector<Input>& inputs, const std::vector<Output>& outputs)
{
    file << "{\"method\":\"" << method << "\",";
    file << "\"function\":\"" << function << "\",";
    WriteArray<Point>(file, "input", inputs);
    WriteArray<Point>(file, "result", outputs);
    file << "}";
}