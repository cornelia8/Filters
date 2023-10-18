#pragma once

#include <string>

#include "components/simple_scene.h"
#include "core/gpu/frame_buffer.h"


namespace m2
{
    class Tema2 : public gfxc::SimpleScene
    {
    public:
        Tema2();
        ~Tema2();

        void Init() override;

    private:
        void FrameStart() override;
        void Update(float deltaTimeSeconds) override;
        void FrameEnd() override;

        void OnInputUpdate(float deltaTime, int mods) override;
        void OnKeyPress(int key, int mods) override;
        void OnKeyRelease(int key, int mods) override;
        void OnMouseMove(int mouseX, int mouseY, int deltaX, int deltaY) override;
        void OnMouseBtnPress(int mouseX, int mouseY, int button, int mods) override;
        void OnMouseBtnRelease(int mouseX, int mouseY, int button, int mods) override;
        void OnMouseScroll(int mouseX, int mouseY, int offsetX, int offsetY) override;
        void OnWindowResize(int width, int height) override;

        void OpenDialog();
        void OnFileSelected(const std::string& fileName);
        void SaveImage(const std::string& fileName);

        // Processing effects
        void img1();
        void img2();
        void GrayScale(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels);
        void Blur(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int blurRadius);
        void Sobel(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels);
        void medianCut(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int nColors);
        void Median(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int blurRadius);
        void nonMaximaSuppression(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int windowSize);
        void hysteresisThresholding(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int low_tresh, int high_tresh);

        int color_distance(float r1, float g1, float b1, float r2, float g2, float b2, int &channel);
        void median_pixel(std::vector<glm::vec3> pixel_list, int start, int end, std::vector<glm::vec3> &color_list, int nColors, int size);

        int BLUR_RADIUS = 3;
        int N_COLORS = 8;
        int WINDOW_SIZE = 5;

    private:
        Texture2D* originalImage;
        Texture2D* processedImage;

        int outputMode;
        bool gpuProcessing;
        bool saveScreenToImage;
    };
}   // namespace m2
