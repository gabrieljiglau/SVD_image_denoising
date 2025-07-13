#pragma once
#include <vector>
#include <optional>
#include <string>


std::optional<std::vector<int>> rgbChannel(std::vector<unsigned char> image, int height, int width, char channel);