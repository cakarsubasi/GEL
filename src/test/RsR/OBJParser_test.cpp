//
// Created by Cem Akarsubasi on 5/21/25.
//

#include <GEL/Util/RawObj.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

using namespace GEL::Util;
using namespace GEL::Util::Combinators;

TEST_CASE("ignore_spaces")
{
    constexpr auto test = "    a";
    constexpr auto test2 = "     ";

    auto view = std::string_view(test);
    auto view2 = std::string_view(test2);
    SUBCASE("simple")
    {
        CHECK_EQ(ignore_spaces(view), true);
        CHECK_EQ(view, "a");
    }
    SUBCASE("ignore all")
    {
        CHECK_EQ(ignore_spaces(view2), true);
        CHECK_EQ(view2, "");
    }
}

TEST_CASE("parse_string")
{
    constexpr auto test = "hello";
    auto view = std::string_view(test);

    SUBCASE("simple")
    {
        CHECK_EQ(parse_string("hel", view), true);
        CHECK_EQ(view, "lo");
    }

}

TEST_CASE("parse_float")
{
    constexpr auto test = "1.5";
    auto view = std::string_view(test);
    SUBCASE("simple")
    {
        CHECK_EQ(parse_float(view), 1.5);
        CHECK_EQ(parse_float(view), std::nullopt);
    }

    constexpr auto test2 = " 1.5";
    auto view2 = std::string_view(test2);
    SUBCASE("leading space")
    {
        CHECK_EQ(parse_float(view2), std::nullopt);
        ignore_spaces(view2);
        CHECK_EQ(parse_float(view2), 1.5);
    }
}

bool float_eq(const double a, const double b, const float eps = 1.0e-6)
{
    return std::abs(a - b) < eps;
}

struct PointCloud {
    std::vector<CGLA::Vec3d> vertices;
    std::vector<CGLA::Vec3d> normals;
};

PointCloud read_obj(std::string& file_path) {
    PointCloud pc;
    std::ifstream file(file_path);
    std::string line;
    if (!file.is_open()) {
        throw std::runtime_error("Error: Unable to open file " + file_path);
    }
    while (std::getline(file, line))
    {
        std::vector<std::string> info;
        size_t pos = 0;
        while ((pos = line.find(' ')) != std::string::npos) {
            info.push_back(line.substr(0, pos));
            line.erase(0, pos + 1);
        }
        info.push_back(line);
        if (info.empty()) {
            continue;
        }
        if (info.at(0) == "v") {
            CGLA::Vec3d vertex(std::stof(info.at(1)),
                std::stof(info.at(2)), std::stof(info.at(3)));
            pc.vertices.emplace_back(vertex);
        }
        if (info.at(0) == "vn") {
            CGLA::Vec3d normal(std::stof(info.at(1)),
                std::stof(info.at(2)), std::stof(info.at(3)));
            pc.normals.emplace_back(normal);
        }
    }
    return pc;
}

TEST_CASE("read raw obj")
{
    auto file_path = "../../../../data/bunny.obj";
    auto file_path_str = std::string(file_path);
    auto robj = read_raw_obj(file_path);
    auto pc = read_obj(file_path_str);
    CHECK_EQ(robj.vertices.size(), pc.vertices.size());
    CHECK_EQ(robj.normals.size(), pc.normals.size());

    std::cout << pc.vertices.size() << std::endl;
}