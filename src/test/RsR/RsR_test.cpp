//
// Created by Cem Akarsubasi on 4/15/25.
//
#include <iostream>
#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/RsR2.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/load.h>

using GEL::HMesh::RSR::point_cloud_to_mesh;

static constexpr auto file_name =
        "../../../../data/PointClouds/Capital_A.obj";
        //    "../../../../data/bunny.obj";
//    "../../../../data/torus.obj";

int test_new() {
    TriMesh input;
    load(file_name, input);
    //assert(success);

    std::cout << "obj vertices: " << input.geometry.no_vertices() << "\n";
    std::cout << "obj normals: " << input.normals.no_vertices() << "\n";

    std::vector<Point> points;
    std::vector<Vector> normals;

    for (auto& point: input.geometry.vertices_copy())
    {
        points.emplace_back(point);
    }
    for (auto& point: input.normals.vertices_copy())
    {
        normals.emplace_back(point);
    }
    normals = {};
    GEL::HMesh::RSR::RsROpts opts;
    opts.k = 30;
    opts.genus = 1;
    opts.r = 20;
    opts.theta = 20;
    HMesh::Manifold output = point_cloud_to_mesh(points, normals, opts);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";
    const bool result = HMesh::obj_save("Bunny_Euc2.obj", output);
    assert(result);
    return 0;
}


int test_old()
{
    TriMesh input;
    load(file_name, input);
    //assert(success);

    std::cout << "obj vertices: " << input.geometry.no_vertices() << "\n";
    std::cout << "obj normals: " << input.normals.no_vertices() << "\n";

    HMesh::Manifold output;
    std::vector<Point> points;
    std::vector<Vector> normals;

    for (auto& point: input.geometry.vertices_copy())
    {
        points.emplace_back(point);
    }
    for (auto& point: input.normals.vertices_copy())
    {
        normals.emplace_back(point);
    }
    normals = {};
    reconstruct_single(output, points, normals, true,  1, 30, 20, 20);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";
    const bool result = HMesh::obj_save("Bunny_Euc2.obj", output);
    assert(result);
    return 0;
}

int main()
{
    test_new();
    return 0;
}