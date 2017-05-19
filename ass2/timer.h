#pragma once
#include <iostream>
#include <ctime>
#include <cstdint>
#include <cctype>
#include <chrono>

/**
 * @brief Timer class for performance measurement
 */
class Timer
{
private:
	std::chrono::high_resolution_clock::time_point start_time;
	int64_t count_down;//in millis
public:
	void start(int64_t millis = 60000){
		count_down = millis;
		reset();
	}
	void reset() {
		start_time =  std::chrono::high_resolution_clock::now();
	}
	// Returns the time in milli-seconds
	float elapsed() {
		using namespace std::chrono;
		return (duration_cast<nanoseconds>(high_resolution_clock::now() - start_time)).count()/1000000.0;
	}
	/**
	 * @brief Timeout
	 * @details Returns true if there is time out
	 */
	bool isTimeout() {
		return count_down <= uint64_t(elapsed());
	}
	/**
	 * @brief Gives present time stamp
	 */
	static int64_t now(){
		using namespace std::chrono;
		return duration_cast<milliseconds>(high_resolution_clock::now().time_since_epoch()).count();
	}
};