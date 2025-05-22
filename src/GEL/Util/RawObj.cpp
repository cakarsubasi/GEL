//
// Created by Cem Akarsubasi on 5/22/25.
//

#include "RawObj.h"

namespace GEL::Util
{

using namespace Combinators;

std::ostream& operator<<(std::ostream& os, const RawObj& obj)
{
    for (auto const& v : obj.vertices) {
        os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }
    for (auto const& vn: obj.normals) {
        os << "vn" << vn[0] << " " << vn[1] << " " << vn[2] << "\n";
    }
    for (auto const& vt: obj.texture_coordinates) {
        os << "vt" << vt[0] << " " << vt[1] << "\n";
    }
    return os;
}

std::istream& operator>>(std::istream& is, RawObj& obj)
{
    obj.vertices.clear();
    std::string line;
    while (auto& more = std::getline(is, line)) {
        std::string_view line_view = line;
        if (auto pos = parse_prefix_then_float_triplet("v", line_view)) {
            obj.vertices.emplace_back(*pos);
        }
        if (auto pos = parse_prefix_then_float_triplet("vn", line_view)) {
            obj.normals.emplace_back(*pos);
        }
        if (auto pos = parse_prefix_then_float_triplet("vt", line_view)) {
            obj.normals.emplace_back(*pos);
        }
    }
    return is;
}

namespace Combinators
{
    bool parse_one(const char c, std::string_view& s) {
    if (s.empty()) {
        return false;
    } else if (s.front() == c) {
        s = s.substr(1);
        return true;
    }
    return false;
}

bool parse_string(const std::string_view& fragment, std::string_view& s) {
    if (s.starts_with(fragment)) {
        s = s.substr(fragment.length());
        return true;
    } else {
        return false;
    }
}

bool ignore_spaces(std::string_view& s) {
    const auto pos = s.find_first_not_of(' ');
    if (pos != std::string::npos) {
        s = s.substr(pos);
    } else {
        s = s.substr(s.length());
    }
    return true;
}

std::optional<double> parse_float(std::string_view& s) {
    double out;
    const auto begin = s.begin();
    const auto end_pos = s.find(' ');
    const auto end = (end_pos == std::string::npos) ? s.end() : s.begin() + end_pos;
    const auto [p, ec] = std::from_chars(begin, end, out);
    if (ec != std::errc()) {
        return std::nullopt;
    } else {
        const auto offset = p - begin;
        s = s.substr(offset);
        return out;
    }
}

std::optional<double> parse_float_ws(std::string_view& s) {
    if (const auto d = parse_float(s)) {
        ignore_spaces(s);
        return d;
    } else {
        return false;
    }
}

std::optional<CGLA::Vec2d> parse_float_doublet_ws(std::string_view& s) {

    auto s_temp = s;
    const auto out1 = parse_float_ws(s_temp);
    if (!out1) return std::nullopt;
    const auto out2 = parse_float_ws(s_temp);
    if (!out2) return std::nullopt;
    s = s_temp;
    return CGLA::Vec2d{*out1, *out2};
}

std::optional<CGLA::Vec3d> parse_float_triplet_ws(std::string_view& s) {

    auto s_temp = s;
    const auto out1 = parse_float_ws(s_temp);
    if (!out1) return std::nullopt;
    const auto out2 = parse_float_ws(s_temp);
    if (!out2) return std::nullopt;
    const auto out3 = parse_float_ws(s_temp);
    if (!out3) return std::nullopt;
    s = s_temp;
    return CGLA::Vec3d{*out1, *out2, *out3};
}

std::optional<CGLA::Vec2d> parse_prefix_then_float_doublet(const std::string_view prefix, std::string_view& s) {
    auto s_temp = s;
    if (!parse_string(prefix, s_temp)) return std::nullopt;
    ignore_spaces(s_temp);
    auto pos = parse_float_doublet_ws(s_temp);
    if (pos) {
        s = s_temp;
        return pos;
    } else {
        return std::nullopt;
    }
}

std::optional<CGLA::Vec3d> parse_prefix_then_float_triplet(const std::string_view prefix, std::string_view& s) {
    auto s_temp = s;
    if (!parse_string(prefix, s_temp)) return std::nullopt;
    ignore_spaces(s_temp);
    auto pos = parse_float_triplet_ws(s_temp);
    if (pos) {
        s = s_temp;
        return pos;
    } else {
        return std::nullopt;
    }
}
}

RawObj read_raw_obj(const std::string& file_path) {
    std::ifstream file(file_path);
    RawObj obj;
    file >> obj;
    return obj;
}

void write_raw_obj(const std::string& file_path, const RawObj& obj)
{
    std::ofstream file(file_path);
    file << obj;
    file.close();
}
}
