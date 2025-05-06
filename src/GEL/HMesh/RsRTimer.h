//
// Created on 4/24/25.
//
#pragma once
#ifndef GEL_HMESH_RSRTIMER_H
#define GEL_HMESH_RSRTIMER_H

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

///
/// TODO: Make this a bit more reliable with a safe API
struct RsR_Timer {
    RsR_Timer() = default;
    std::vector<long> times;
    std::vector<std::string> descriptions;
    std::vector<std::chrono::high_resolution_clock::time_point> starts;
    std::vector<std::chrono::high_resolution_clock::time_point> ends;

    void create(const std::string& name) {
        times.push_back(0);
        descriptions.push_back(name);
        //idx_map.insert(std::pair<std::string, int>(name, times.size() - 1));
        starts.push_back(std::chrono::high_resolution_clock::now());
        ends.push_back(std::chrono::high_resolution_clock::now());
    }

    void start(const std::string& name) {
        //int idx = idx_map[name];
        int idx = -1;
        for (int i = 0; i < descriptions.size(); i++) {
            if (descriptions[i] == name) {
                idx = i;
                break;
            }
        }
        starts[idx] = std::chrono::high_resolution_clock::now();
    }

    void end(const std::string& name) {
        //int idx = idx_map[name];
        int idx = -1;
        for (int i = 0; i < descriptions.size(); i++) {
            if (descriptions[i] == name) {
                idx = i;
                break;
            }
        }
        ends[idx] = std::chrono::high_resolution_clock::now();
        times[idx] +=
            std::chrono::duration_cast<std::chrono::milliseconds>(ends[idx] - starts[idx]).count();
    }

    void show() const
    {
        std::cout << "Time Statistics" << std::endl;
        std::cout << std::string(20, '=') << std::endl;
        for (int i = 0; i < times.size(); i++) {
            std::cout << "Spent " << times[i]
                << " milliseconds on " << descriptions[i] << std::endl;
        }
    }

    long log(const std::string& name) const
    {
        int idx = -1;
        for (int i = 0; i < descriptions.size(); i++) {
            if (descriptions[i] == name) {
                idx = i;
                break;
            }
        }
        return times[idx];
    }
};

#endif //GEL_HMESH_RSRTIMER_H
