//
// Created by Cem Akarsubasi on 5/22/25.
//
#ifndef RAWOBJ_H
#define RAWOBJ_H

#include <iostream>
#include <fstream>
#include <charconv>
#include <vector>
#include <optional>

#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Vec2d.h>

namespace GEL::Util
{

/// Parser combinators for easy and low-vulnerability parsing
///
namespace Combinators
{
    bool parse_one(char c, std::string_view& s);

    bool parse_string(const std::string_view& fragment, std::string_view& s);

    bool ignore_spaces(std::string_view& s);

    std::optional<double> parse_float(std::string_view& s);

    std::optional<double> parse_float_ws(std::string_view& s);

    std::optional<CGLA::Vec2d> parse_float_doublet_ws(std::string_view& s);

    std::optional<CGLA::Vec3d> parse_float_triplet_ws(std::string_view& s);

    std::optional<CGLA::Vec2d> parse_prefix_then_float_doublet(std::string_view prefix, std::string_view& s);

    std::optional<CGLA::Vec3d> parse_prefix_then_float_triplet(std::string_view prefix, std::string_view& s);
}

/// Raw Wavefront Obj type for
///
struct RawObj {
    std::vector<CGLA::Vec3d> vertices;
    std::vector<CGLA::Vec3d> normals;
    std::vector<CGLA::Vec2d> texture_coordinates;
    // TODO:
    // struct Face {
    //   vertex_id;
    //   tcoord_id;
    //   normal_id;
    // }
    // std::vector<Face> faces;

    friend std::ostream& operator<<(std::ostream& os, const RawObj& obj);
    friend std::istream& operator>>(std::istream& is, RawObj& obj);
};

RawObj read_raw_obj(const std::string& file_path);

void write_raw_obj(const std::string& file_path, const RawObj& obj);
}

#endif //RAWOBJ_H
