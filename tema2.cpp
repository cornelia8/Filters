#include "lab_m2/tema2/tema2.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>

#include "pfd/portable-file-dialogs.h"

using namespace std;
using namespace m2;


/*
 *  To find out more about `FrameStart`, `Update`, `FrameEnd`
 *  and the order in which they are called, see `world.cpp`.
 */

Tema2::Tema2()
{
    outputMode = 0;
    gpuProcessing = false;
    saveScreenToImage = false;
    window->SetSize(600, 600);
}

Tema2::~Tema2()
{
}


void Tema2::Init()
{
    // Load default texture fore imagine processing
    originalImage = TextureManager::LoadTexture(PATH_JOIN(window->props.selfDir, RESOURCE_PATH::TEXTURES, "cube", "pos_x.png"), nullptr, "image", true, true);
    processedImage = TextureManager::LoadTexture(PATH_JOIN(window->props.selfDir, RESOURCE_PATH::TEXTURES, "cube", "pos_x.png"), nullptr, "newImage", true, true);

    {
        Mesh* mesh = new Mesh("quad");
        mesh->LoadMesh(PATH_JOIN(window->props.selfDir, RESOURCE_PATH::MODELS, "primitives"), "quad.obj");
        mesh->UseMaterials(false);
        meshes[mesh->GetMeshID()] = mesh;
    }

    std::string shaderPath = PATH_JOIN(window->props.selfDir, SOURCE_PATH::M2, "tema2", "shaders");

    // Create a shader program for particle system
    {
        Shader* shader = new Shader("ImageProcessing");
        shader->AddShader(PATH_JOIN(shaderPath, "VertexShader.glsl"), GL_VERTEX_SHADER);
        shader->AddShader(PATH_JOIN(shaderPath, "FragmentShader.glsl"), GL_FRAGMENT_SHADER);

        shader->CreateAndLink();
        shaders[shader->GetName()] = shader;
    }
}


void Tema2::FrameStart()
{
}

void Tema2::Update(float deltaTimeSeconds)
{
    ClearScreen();

    auto shader = shaders["ImageProcessing"];
    shader->Use();

    if (saveScreenToImage)
    {
        window->SetSize(originalImage->GetWidth(), originalImage->GetHeight());
    }

    int flip_loc = shader->GetUniformLocation("flipVertical");
    glUniform1i(flip_loc, saveScreenToImage ? 0 : 1);

    int screenSize_loc = shader->GetUniformLocation("screenSize");
    glm::ivec2 resolution = window->GetResolution();
    glUniform2i(screenSize_loc, resolution.x, resolution.y);

    int outputMode_loc = shader->GetUniformLocation("outputMode");
    glUniform1i(outputMode_loc, outputMode);

    int locTexture = shader->GetUniformLocation("textureImage");
    glUniform1i(locTexture, 0);

    auto textureImage = (gpuProcessing == true) ? originalImage : processedImage;
    textureImage->BindToTextureUnit(GL_TEXTURE0);

    RenderMesh(meshes["quad"], shader, glm::mat4(1));

    if (saveScreenToImage)
    {
        saveScreenToImage = false;

        GLenum format = GL_RGB;
        if (originalImage->GetNrChannels() == 4)
        {
            format = GL_RGBA;
        }

        glReadPixels(0, 0, originalImage->GetWidth(), originalImage->GetHeight(), format, GL_UNSIGNED_BYTE, processedImage->GetImageData());
        processedImage->UploadNewData(processedImage->GetImageData());
        SaveImage("shader_processing_" + std::to_string(outputMode));

        float aspectRatio = static_cast<float>(originalImage->GetWidth()) / originalImage->GetHeight();
        window->SetSize(static_cast<int>(600 * aspectRatio), 600);
    }
}

void Tema2::FrameEnd()
{
    DrawCoordinateSystem();
}

void Tema2::OpenDialog()
{
    std::vector<std::string> filters =
    {
        "Image Files", "*.png *.jpg *.jpeg *.bmp",
        "All Files", "*"
    };

    auto selection = pfd::open_file("Select a file", ".", filters).result();
    if (!selection.empty())
    {
        std::cout << "User selected file " << selection[0] << "\n";
        OnFileSelected(selection[0]);
    }
}

void Tema2::OnFileSelected(const std::string& fileName)
{
    if (fileName.size())
    {
        std::cout << fileName << endl;
        originalImage = TextureManager::LoadTexture(fileName, nullptr, "image", true, true);
        processedImage = TextureManager::LoadTexture(fileName, nullptr, "newImage", true, true);

        float aspectRatio = static_cast<float>(originalImage->GetWidth()) / originalImage->GetHeight();
        window->SetSize(static_cast<int>(600 * aspectRatio), 600);
    }
}

void Tema2::SaveImage(const std::string& fileName)
{
    cout << "Saving image! ";
    processedImage->SaveToFile((fileName + ".png").c_str());
    cout << "[Done]" << endl;
}

void Tema2::img1()
{
    unsigned int channels = originalImage->GetNrChannels();
    unsigned char* data = originalImage->GetImageData();
    unsigned char* newData = processedImage->GetImageData();

    if (channels < 3)
        return;

    glm::ivec2 imageSize = glm::ivec2(originalImage->GetWidth(), originalImage->GetHeight());

    unsigned char* data1 = (unsigned char*)calloc(channels * imageSize.y * imageSize.x, sizeof(unsigned char));
    unsigned char* data2 = (unsigned char*)calloc(channels * imageSize.y * imageSize.x, sizeof(unsigned char));

    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            int offset = channels * (i * imageSize.x + j);
            // copy original
            data1[offset] = data[offset];
            data1[offset + 1] = data[offset + 1];
            data1[offset + 2] = data[offset + 2];
        }
    }

    // grayscale (in = data1, out = data2)
    GrayScale(data1, data2, imageSize, channels);
    // blur (in = data2, out = data1)
    Blur(data2, data1, imageSize, channels, BLUR_RADIUS);
    // sobel (in = data1, out = data2)
    Sobel(data1, data2, imageSize, channels);
    // non-maxima supression (in = data2, out = data1)
    nonMaximaSuppression(data2, data1, imageSize, channels, WINDOW_SIZE);
    // hysteresis thresholding (in = data1, out = data2)

    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            int offset = channels * (i * imageSize.x + j);
            // copy result
            newData[offset] = data2[offset];
            newData[offset + 1] = data2[offset + 1];
            newData[offset + 2] = data2[offset + 2];
        }
    }

    free(data1);
    free(data2);

    processedImage->UploadNewData(newData);
}

void Tema2::img2()
{
    unsigned int channels = originalImage->GetNrChannels();
    unsigned char* data = originalImage->GetImageData();
    unsigned char* newData = processedImage->GetImageData();

    if (channels < 3)
        return;

    glm::ivec2 imageSize = glm::ivec2(originalImage->GetWidth(), originalImage->GetHeight());

    unsigned char* data1 = (unsigned char*)calloc(channels * imageSize.y * imageSize.x, sizeof(unsigned char));
    unsigned char* data2 = (unsigned char*)calloc(channels * imageSize.y * imageSize.x, sizeof(unsigned char));

    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            int offset = channels * (i * imageSize.x + j);
            // copy original
            data1[offset] = data[offset];
            data1[offset + 1] = data[offset + 1];
            data1[offset + 2] = data[offset + 2];
        }
    }

    // median cut (in = data1, out = data2)
    medianCut(data1, data2, imageSize, channels, N_COLORS);
    // median filter (in = data2, out = data1)
    Median(data2, data1, imageSize, channels, BLUR_RADIUS);

    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            int offset = channels * (i * imageSize.x + j);
            // copy result
            newData[offset] = data1[offset];
            newData[offset + 1] = data1[offset + 1];
            newData[offset + 2] = data1[offset + 2];
        }
    }

    free(data1);
    free(data2);

    processedImage->UploadNewData(newData);
}

void Tema2::GrayScale(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels) 
{
    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            int offset = channels * (i * imageSize.x + j);

            // Reset save image data
            char value = static_cast<char>(data[offset + 0] * 0.2f + data[offset + 1] * 0.71f + data[offset + 2] * 0.07);
            memset(&newData[offset], value, 3);
        }
    }
}

void Tema2::Blur(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int blurRadius)
{
    int offset, offset_neighb;

    for (int i = blurRadius; i < imageSize.y - blurRadius; i++)
    {
        for (int j = blurRadius; j < imageSize.x - blurRadius; j++)
        {
            int offset = channels * (i * imageSize.x + j);
            int sum_r = 0;
            int sum_g = 0;
            int sum_b = 0;
            for (int y = i - blurRadius; y <= i + blurRadius; y++) {
                for (int x = j - blurRadius; x <= j + blurRadius; x++) {
                    offset_neighb = channels * (y * imageSize.x + x);

                    sum_r += data[offset_neighb + 0];
                    sum_g += data[offset_neighb + 1];
                    sum_b += data[offset_neighb + 2];
                }
            }
            // Reset save image data

            newData[offset + 0] = sum_r / pow((2 * blurRadius + 1), 2);
            newData[offset + 1] = sum_g / pow((2 * blurRadius + 1), 2);
            newData[offset + 2] = sum_b / pow((2 * blurRadius + 1), 2);
        }
    }
}

void Tema2::Sobel(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels)
{
    int offset, offset_neighb;
    int mask_x[] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
    int mask_y[] = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };
    int mask_idx = 0;
    int sum_x = 0;
    int sum_y = 0;

    for (int y = 1; y < imageSize.y - 1; y++)
    {
        for (int x = 1; x < imageSize.x - 1; x++)
        {
            offset = (y * imageSize.x + x) * channels;
            mask_idx = 0;
            sum_x = 0;
            sum_y = 0;

            for (int j = y - 1; j <= y + 1; j++)
            {
                for (int i = x - 1; i <= x + 1; i++)
                {
                    offset_neighb = (j * imageSize.x + i) * 3;

                    sum_x += data[offset_neighb] * mask_x[mask_idx];
                    sum_y += data[offset_neighb] * mask_y[mask_idx];
                    mask_idx++;
                }
            }
            newData[offset] = sqrt(sum_x * sum_x + sum_y * sum_y);
            newData[offset + 1] = sqrt(sum_x * sum_x + sum_y * sum_y);
            newData[offset + 2] = sqrt(sum_x * sum_x + sum_y * sum_y);
        }
    }
}

int Tema2::color_distance(float r1, float g1, float b1, float r2, float g2, float b2, int &channel)
{
    int dif = abs(r1 - r2);
    channel = 1;
    if (abs(g1 - g2) > dif)
    {
        dif = abs(g1 - g2);
        channel = 2;
    }
    if (abs(b1 - b2) > dif)
    {
        dif = abs(b1 - b2);
        channel = 2;
    }

    return (int) (abs(r1 - r2) + abs(g1 - g2) + abs(b1 - b2));
}

void Tema2::median_pixel(std::vector<glm::vec3> pixel_list, int start, int end, std::vector<glm::vec3> &color_list, int nColors, int size)
{
    glm::vec3 median_pixel1 = glm::vec3(0, 0, 0);
    glm::vec3 median_pixel2 = glm::vec3(0, 0, 0);

    // color_distance(pixel1.x, pixel1.y, pixel1. z, pixel2.x, pixel2.y, pixel2.z, channel);
    int max_dif = 0;
    int channel = 0;
    glm::vec3 pixel1, pixel2;
    glm::vec3 min_pixel = pixel_list.at(start), max_pixel = pixel_list.at(start);
    int min_sum, max_sum, index_min = start, index_max = start;
    min_sum = min_pixel.x + min_pixel.y + min_pixel.z;
    max_sum = min_sum;

    // cautarea diferentei maxime si canalului cel mai proeminent
    for (int i = start+1; i < end - 1; i++)
    {
        pixel1 = pixel_list.at(i);
        int sum = pixel1.x + pixel1.y + pixel1.z;
        if (sum < min_sum)
        {
            min_sum = sum;
            index_min = i;
        }
        if (sum > max_sum)
        {
            max_sum = sum;
            index_max = i;
        }
    }
    pixel1 = pixel_list.at(index_min);
    pixel2 = pixel_list.at(index_max);
    max_dif = color_distance(pixel1.x, pixel1.y, pixel1.z, pixel2.x, pixel2.y, pixel2.z, channel);

    // printf("sorting after channel %d\n", channel);

    // sortarea listei crescator dupa canalul ales
    switch(channel)
    {
    case 1:
        std::sort(std::next(pixel_list.begin(), start), std::next(pixel_list.begin(), end), 
            [](const glm::vec3& a, const glm::vec3& b) -> bool { return a.x < b.x; });
        break;
    case 2:
        std::sort(std::next(pixel_list.begin(), start), std::next(pixel_list.begin(), end),
            [](const glm::vec3& a, const glm::vec3& b) -> bool { return a.y < b.y; });
        break;
    case 3:
        std::sort(std::next(pixel_list.begin(), start), std::next(pixel_list.begin(), end),
            [](const glm::vec3& a, const glm::vec3& b) -> bool { return a.z < b.z; });
        break;
    default:
        std::sort(std::next(pixel_list.begin(), start), std::next(pixel_list.begin(), end),
            [](const glm::vec3& a, const glm::vec3& b) -> bool { return a.x < b.x; });
        break;
    }

    // impartirea pe 2 intervale si calculul pixelului mediu pe fiecare interval
    int middle = (end + start) / 2;
    int nr = (end - start) / 2;

    for (int i = start; i < middle; i++)
    {
        median_pixel1.x = median_pixel1.x + pixel_list.at(i).x;
        median_pixel1.y = median_pixel1.y + pixel_list.at(i).y;
        median_pixel1.z = median_pixel1.z + pixel_list.at(i).z;
    }

    median_pixel1.x = median_pixel1.x / nr;
    median_pixel1.y = median_pixel1.y / nr;
    median_pixel1.z = median_pixel1.z / nr;

    for (int i = middle; i < end; i++)
    {
        median_pixel2.x = median_pixel2.x + pixel_list.at(i).x;
        median_pixel2.y = median_pixel2.y + pixel_list.at(i).y;
        median_pixel2.z = median_pixel2.z + pixel_list.at(i).z;
    }

    median_pixel2.x = median_pixel2.x / nr;
    median_pixel2.y = median_pixel2.y / nr;
    median_pixel2.z = median_pixel2.z / nr;

    if (size*2 < nColors)
    {
        median_pixel(pixel_list, start, middle, color_list, nColors, size*2);
        median_pixel(pixel_list, middle, end, color_list, nColors, size*2);
    }

    if (size*2 == nColors)
    {
        color_list.push_back(median_pixel1);
        color_list.push_back(median_pixel2);
    }
}

void Tema2::medianCut(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int nColors)
{
    std::vector<glm::vec3> pixel_list;
    int offset;

    // lista de vectori din imaginea originala
    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            offset = channels * (i * imageSize.x + j);
            pixel_list.push_back(glm::vec3(data[offset], data[offset+1], data[offset+2]));
        }
    }

    // popularea listei de culori
    int size = 1;
    std::vector<glm::vec3> color_list;
    median_pixel(pixel_list, 0, pixel_list.size() / 2, color_list, nColors, size*2);
    median_pixel(pixel_list, pixel_list.size() / 2, pixel_list.size(), color_list, nColors, size*2);

    // printf("Size of color_list is %d\n", color_list.size());

    // inlocuirea culorilor
    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            offset = channels * (i * imageSize.x + j);
            glm::vec3 old_pixel = glm::vec3(data[offset], data[offset + 1], data[offset + 2]);

            // cautarea pixelului cel mai apropiat
            int min_dif = -1;
            int index = -1;
            for (int k = 0; k < color_list.size(); k++)
            {
                int temp_channel = 0;
                int dif = color_distance(old_pixel.x, old_pixel.y, old_pixel.z, color_list.at(k).x, color_list.at(k).y, color_list.at(k).z, temp_channel);

                if ((min_dif == -1) || (dif < min_dif))
                {
                    min_dif = dif;
                    index = k;
                }
            }

            newData[offset] = (int) color_list.at(index).x;
            newData[offset + 1] = (int) color_list.at(index).y;
            newData[offset + 2] = (int) color_list.at(index).z;
        }
    }
}

void Tema2::Median(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int blurRadius)
{
    int offset;

    for (int i = blurRadius; i < imageSize.y - blurRadius; i++)
    {
        for (int j = blurRadius; j < imageSize.x - blurRadius; j++)
        {
            offset = channels * (i * imageSize.x + j);
            vector<glm::vec3> v;
            int offset_neighb;

            for (int x = i - blurRadius; x <= i + blurRadius; x++)
            {
                for (int y = j - blurRadius; y <= j + blurRadius; y++)
                {
                    offset_neighb = channels * (x * imageSize.x + y);

                    v.push_back(glm::vec3(data[offset_neighb], data[offset_neighb + 1], data[offset_neighb + 2]));
                }
            }

            std::sort(v.begin(), v.end(), [](const glm::vec3& a, const glm::vec3& b) -> bool 
                { return (a.x * 0.2f + a.y * 0.71f + a.z * 0.07f) < (b.x * 0.2f + b.y * 0.71f + b.z * 0.07f); });
            glm::vec3 pixel = v.at(v.size() / 2);

            newData[offset] = pixel.x;
            newData[offset + 1] = pixel.y;
            newData[offset + 2] = pixel.z;
        }
    }

    for (int i = 0; i < blurRadius; i++)
    {
        for (int j = 0; j < blurRadius; j++)
        {
            offset = channels * (i * imageSize.x + j);

            newData[offset] = data[offset];
            newData[offset + 1] = data[offset + 1];
            newData[offset + 2] = data[offset + 2];
        }
    }

    for (int i = imageSize.y - blurRadius; i < imageSize.y; i++)
    {
        for (int j = imageSize.x - blurRadius; j < imageSize.x; j++)
        {
            offset = channels * (i * imageSize.x + j);

            newData[offset] = data[offset];
            newData[offset + 1] = data[offset + 1];
            newData[offset + 2] = data[offset + 2];
        }
    }
}

void Tema2::nonMaximaSuppression(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int windowSize)
{
    int offset;

    for (int i = 0; i < imageSize.y; i++)
    {
        for (int j = 0; j < imageSize.x; j++)
        {
            offset = channels * (i * imageSize.x + j);
            int max_val = 0;
            int max_r = 0, max_g = 0, max_b = 0;
            int offset_neighb;

            for (int y = i - windowSize; y <= i + windowSize; y++)
            {
                for (int x = j -windowSize; x <= j + windowSize; x++)
                {
                    if (y < 0 || x < 0 || y > imageSize.y || x > imageSize.x)
                        continue;

                    offset_neighb = channels * (x * imageSize.x + y);
                    int r, g, b;
                    r = data[offset_neighb];
                    g = data[offset_neighb + 1];
                    b = data[offset_neighb + 2];

                    if (r + g + b > max_val)
                    {
                        max_val = r + g + b;
                        max_r = r;
                        max_g = g;
                        max_b = b;
                    }
                }
            }

            if (data[offset] == max_r && data[offset + 1] == max_g && data[offset + 2] == max_b)
            {
                newData[offset] = data[offset];
                newData[offset + 1] = data[offset + 1];
                newData[offset + 2] = data[offset + 2];
            }
            else 
            {
                newData[offset] = 0;
                newData[offset + 1] = 0;
                newData[offset + 2] = 0;
            }
        }
    }
}

void Tema2::hysteresisThresholding(unsigned char* data, unsigned char* newData, glm::ivec2 imageSize, unsigned int channels, int low_tresh, int high_tresh)
{
    int offset;

    for (int i = 1; i < imageSize.y - 1; i++)
    {
        for (int j = 1; j < imageSize.x - 1; j++)
        {
            offset = channels * (i * imageSize.x + j);
        }
    }
}

/*
 *  These are callback functions. To find more about callbacks and
 *  how they behave, see `input_controller.h`.
 */

void Tema2::OnInputUpdate(float deltaTime, int mods)
{
    // Treat continuous update based on input
}

void Tema2::OnKeyPress(int key, int mods)
{
    // Add key press event
    if (key == GLFW_KEY_F || key == GLFW_KEY_ENTER || key == GLFW_KEY_SPACE)
    {
        OpenDialog();
    }

    if (key == GLFW_KEY_E)
    {
        gpuProcessing = !gpuProcessing;
        if (gpuProcessing == false)
        {
            outputMode = 0;
        }
        cout << "Processing on GPU: " << (gpuProcessing ? "true" : "false") << endl;
    }

    if (key - GLFW_KEY_0 >= 0 && key < GLFW_KEY_7)
    {
        outputMode = key - GLFW_KEY_0;

        if (gpuProcessing == false)
        {
            outputMode = key - GLFW_KEY_0;

            if (gpuProcessing == false)
            {
                switch (outputMode)
                {
                case 1:
                    img1();
                    break;
                case 2:
                    img2();
                    break;
                default:
                    processedImage->UploadNewData(originalImage->GetImageData());
                    break;
                }
            }
        }
    }

    if (key == GLFW_KEY_S && mods & GLFW_MOD_CONTROL)
    {
        if (!gpuProcessing)
        {
            SaveImage("processCPU_" + std::to_string(outputMode));
        }
        else {
            saveScreenToImage = true;
        }
    }
}


void Tema2::OnKeyRelease(int key, int mods)
{
    // Add key release event
}


void Tema2::OnMouseMove(int mouseX, int mouseY, int deltaX, int deltaY)
{
    // Add mouse move event
}


void Tema2::OnMouseBtnPress(int mouseX, int mouseY, int button, int mods)
{
    // Add mouse button press event
}


void Tema2::OnMouseBtnRelease(int mouseX, int mouseY, int button, int mods)
{
    // Add mouse button release event
}


void Tema2::OnMouseScroll(int mouseX, int mouseY, int offsetX, int offsetY)
{
    // Treat mouse scroll event
}


void Tema2::OnWindowResize(int width, int height)
{
    // Treat window resize event
}
