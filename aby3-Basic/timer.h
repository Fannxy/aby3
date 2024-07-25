#pragma once
#include <chrono>
#include <iostream>
#include <unordered_map>
#include <map>
#include <string>
#include "debug.h"

class Timer{
public:

    // Singleton pattern.
    static Timer& getInstance() {
        static Timer instance; // Guaranteed to be destroyed.
                               // Instantiated on first use.
        return instance;
    }

    Timer(Timer const&) = delete; // you can not copy Timer instance.
    void operator=(Timer const&) = delete; // you can not assign Timer instance.

    void start(const std::string& key) {
        timestamps[key] = std::chrono::high_resolution_clock::now();
    }

    void end(const std::string& key) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto start_time = timestamps[key];
        durations[key] = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }

    std::string get_key(const std::string& key) {
        int count = keycounts[key]++;
        return key + "-" + std::to_string(count);
    }

    void print(const std::string& key, const std::string& unit = "microseconds", std::ostream& os = std::cout) {
        // check if the key exists
        if(durations.find(key) == durations.end()){
            THROW_RUNTIME_ERROR("Key: " + key + " not found in the timer");
        }
        
        double duration = durations[key];
        if (unit == "milliseconds") {
            duration /= 1000;
        } else if (unit == "seconds") {
            duration /= 1000000;
        } else if (unit == "minutes") {
            duration /= 60000000;
        }
        os << "Time taken by " << key << ": " << duration << " " << unit << std::endl;
    }

    void print(const std::string& unit = "microseconds", std::ostream& os = std::cout) {
        for(auto& kv : durations){
            print(kv.first, unit, os);
        }
    }

    void print_from_total(const std::string& key, const std::string& unit = "microseconds", std::ostream& os = std::cout){
        // check if the key exists
        if(totalTimes.find(key) == totalTimes.end()){
            THROW_RUNTIME_ERROR("Key: " + key + " not found in the timer");
        }
        
        double duration = totalTimes[key];
        if (unit == "milliseconds") {
            duration /= 1000;
        } else if (unit == "seconds") {
            duration /= 1000000;
        } else if (unit == "minutes") {
            duration /= 60000000;
        }
        os << "Time taken by " << key << ": " << duration << " " << unit << std::endl;
    }

    void print_total(const std::string& unit = "microseconds", std::ostream& os = std::cout){
        accumulate();
        for(auto& kv : totalTimes){
            print_from_total(kv.first, unit, os);
        }
    }

    void print_records(const std::string& key_prefix, const std::string& unit = "microseconds", std::ostream& os = std::cout){
        for(auto& kv : durations){
            if(kv.first.substr(0, key_prefix.size()) == key_prefix){
                print(kv.first, unit, os);
            }
        }
    }

    void accumulate() {
        for(auto& kv : durations){
            std::string time_prefix = get_prefix(kv.first);
            if(totalTimes.find(time_prefix) == totalTimes.end()){
                totalTimes[time_prefix] = kv.second;
            }
            else{
                totalTimes[time_prefix] += kv.second;
            }
        }
        return;
    }

    void clear_records(){
        timestamps.clear();
        durations.clear();
        keycounts.clear();
        totalTimes.clear();
    }

private:
    Timer() {} // you can not construct Timer instance from outside.
    std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> timestamps;
    std::map<std::string, double> durations;
    std::unordered_map<std::string, int> keycounts;
    std::unordered_map<std::string, double> totalTimes;

    std::string get_prefix(const std::string& key) {
        size_t pos = key.find_last_of('-');
        if (pos != std::string::npos) {
            return key.substr(0, pos);
        } else {
            return key;
        }
    }
};

class CommunicationMeter {
public:
    // Singleton pattern.
    static CommunicationMeter& getInstance() {
        static CommunicationMeter instance; // Guaranteed to be destroyed.
                                            // Instantiated on first use.
        return instance;
    }

    CommunicationMeter(CommunicationMeter const&) = delete; // you can not copy CommunicationMeter instance.
    void operator=(CommunicationMeter const&) = delete; // you can not assign CommunicationMeter instance.

    void start(const std::string& key, uint64_t start_communication) {
        start_communications[key] = start_communication;
    }

    void end(const std::string& key, uint64_t end_communication) {
        int start_communication = start_communications[key];
        communications[key] = end_communication - start_communication;
    }

    int getCommunication(const std::string& key) {
        return communications[key];
    }

    void print(const std::string& key, const std::string& unit = "MB", std::ostream& os = std::cout) {
        if(communications.find(key) == communications.end()){
            THROW_RUNTIME_ERROR("Key: " + key + " not found in the timer");
        }
        
        double comm = communications[key];
        if (unit == "KB") {
            comm /= 1000;
        } else if (unit == "MB") {
            comm /= 1000000;
        } else if (unit == "GB") {
            comm /= 1000000000;
        }
        os << "Communicaitions of by " << key << ": " << comm << " " << unit << std::endl;
    }

    void print_records(const std::string& key_prefix, const std::string& unit = "MB", std::ostream& os = std::cout){
        for(auto& kv : communications){
            if(kv.first.substr(0, key_prefix.size()) == key_prefix){
                print(kv.first, unit, os);
            }
        }
    }

    void print_total(const std::string& unit = "MB", std::ostream& os = std::cout){
        for(auto& kv : communications){
            print(kv.first, unit, os);
        }
    }

    void print_total_per_party(const std::string& unit = "MB", std::ostream& os = std::cout){
        for(auto& kv : totalCommunications){
            std::string key = kv.first;
            double total_comm = 0;
            uint64_t comm1, comm2, comm3;
            std::tie(comm1, comm2, comm3) = kv.second; 
            std::array<uint64_t, 3> comms = {comm1, comm2, comm3};
            for(size_t i=0; i<3; i++){
                double party_comm = comms[i];
                if (unit == "KB") {
                    party_comm /= 1000;
                } else if (unit == "MB") {
                    party_comm /= 1000000;
                } else if (unit == "GB") {
                    party_comm /= 1000000000;
                }
                total_comm += party_comm;
                os << "Communications of " << key << " Party - " << std::to_string(i) << " : " << party_comm << " " << unit << std::endl;
            }
            os << "Total Communications of " << key << " : " << total_comm << " " << unit << std::endl;
        }
        return;
    }

    CommunicationMeter() {}

    std::map<std::string, uint64_t> start_communications;
    std::map<std::string, uint64_t> communications;
    std::map<std::string, std::tuple<uint64_t, uint64_t, uint64_t>> totalCommunications;
};