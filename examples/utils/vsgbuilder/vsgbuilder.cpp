#include <vsg/all.h>

#ifdef vsgXchange_FOUND
#    include <vsgXchange/all.h>
#endif

#include <iostream>

struct FindMatchingBufferInfo : public vsg::Visitor
{
    vsg::ref_ptr<vsg::Data> data;
    std::set<vsg::ref_ptr<vsg::BufferInfo>> matches;

    FindMatchingBufferInfo(vsg::ref_ptr<vsg::Data> in_data) : data(in_data) {}


    void apply(vsg::Object& object) override
    {
        object.traverse(*this);
    }

    void apply(vsg::VertexIndexDraw& vid) override
    {
        for(auto& bufferInfo : vid.arrays)
        {
            if (bufferInfo->data == data) matches.insert(bufferInfo);
        }
    }
};


int main(int argc, char** argv)
{
    // set up defaults and read command line arguments to override them
    auto options = vsg::Options::create();
    options->paths = vsg::getEnvPaths("VSG_FILE_PATH");
    options->objectCache = vsg::ObjectCache::create();

    auto windowTraits = vsg::WindowTraits::create();
    windowTraits->windowTitle = "vsgbuilder";

    auto builder = vsg::Builder::create();
    builder->options = options;

    // set up defaults and read command line arguments to override them
    vsg::CommandLine arguments(&argc, argv);
    windowTraits->debugLayer = arguments.read({"--debug", "-d"});
    windowTraits->apiDumpLayer = arguments.read({"--api", "-a"});

    arguments.read("--screen", windowTraits->screenNum);
    arguments.read("--display", windowTraits->display);
    auto numFrames = arguments.value(-1, "-f");
    if (arguments.read({"--fullscreen", "--fs"})) windowTraits->fullscreen = true;
    if (arguments.read({"--window", "-w"}, windowTraits->width, windowTraits->height)) { windowTraits->fullscreen = false; }
    if (arguments.read("--IMMEDIATE")) windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_IMMEDIATE_KHR;
    if (arguments.read("--double-buffer")) windowTraits->swapchainPreferences.imageCount = 2;
    if (arguments.read("--triple-buffer")) windowTraits->swapchainPreferences.imageCount = 3; // defaul
    if (arguments.read("-t"))
    {
        windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_IMMEDIATE_KHR;
        windowTraits->width = 192, windowTraits->height = 108;
        windowTraits->decoration = false;
    }

    bool floatColors = !arguments.read("--ubvec4-colors");

    bool box = arguments.read("--box");
    bool capsule = arguments.read("--capsule");
    bool cone = arguments.read("--cone");
    bool cylinder = arguments.read("--cylinder");
    bool quad = arguments.read("--quad");
    bool sphere = arguments.read("--sphere");

    bool dynamicUpdates = arguments.read("--dynamic");
    size_t frameToWait = arguments.value<size_t>(2, {"--frames-to-wait", "--f2w"});
    auto waitTimeout = arguments.value<uint64_t>(50000000, {"--timeout", "--to"});

    if (!(box || sphere || cone || capsule || quad || cylinder))
    {
        box = true;
        capsule  = true;
        cone = true;
        cylinder = true;
        quad = true;
        sphere = true;
    }

    auto numVertices = arguments.value<uint32_t>(0, "-n");

    vsg::Path textureFile = arguments.value(vsg::Path{}, {"-i", "--image"});

    if (arguments.errors()) return arguments.writeErrorMessages(std::cerr);

#ifdef vsgXchange_all
    // add vsgXchange's support for reading and writing 3rd party file formats
    options->add(vsgXchange::all::create());
#endif

    auto scene = vsg::Group::create();

    vsg::dvec3 centre = {0.0, 0.0, 0.0};
    double radius = 1.0;

    vsg::GeometryInfo info;
    {
        info.dx.set(1.0f, 0.0f, 0.0f);
        info.dy.set(0.0f, 1.0f, 0.0f);
        info.dz.set(0.0f, 0.0f, 1.0f);

        if (!textureFile.empty()) info.image = vsg::read_cast<vsg::Data>(textureFile, options);

        vsg::dbox bound;

        if (numVertices>0)
        {
            float w = std::pow(float(numVertices), 0.33f) * 2.0f;
            info.positions = vsg::vec3Array::create(numVertices);
            for(auto& v : *(info.positions))
            {
                v.set(w * (float(std::rand()) / float(RAND_MAX) - 0.5f),
                      w * (float(std::rand()) / float(RAND_MAX) - 0.5f),
                      w * (float(std::rand()) / float(RAND_MAX) - 0.5f));
            }

            radius += (0.5 * sqrt(3.0) * w) ;

            if (floatColors)
            {
                auto colors = vsg::vec4Array::create(info.positions->size());
                info.colors = colors;
                for(auto& c : *(colors))
                {
                    c.set(float(std::rand()) / float(RAND_MAX), float(std::rand()) / float(RAND_MAX), float(std::rand()) / float(RAND_MAX), 1.0f);
                }
            }
            else
            {
                auto colors = vsg::ubvec4Array::create(info.positions->size());
                info.colors = colors;
                for(auto& c : *(colors))
                {
                    c.set(uint8_t(255.0 * float(std::rand()) / float(RAND_MAX)), uint8_t(255.0 * float(std::rand()) / float(RAND_MAX)), uint8_t(255.0 * float(std::rand()) / float(RAND_MAX)), 255);
                }
            }

        }

        if (box)
        {
            scene->addChild(builder->createBox(info));
            bound.add(info.position);
            info.position += info.dx * 1.5f;
        }

        if (sphere)
        {
            scene->addChild(builder->createSphere(info));
            bound.add(info.position);
            info.position += info.dx * 1.5f;
        }

        if (quad)
        {
            scene->addChild(builder->createQuad(info));
            bound.add(info.position);
            info.position += info.dx * 1.5f;
        }

        if (cylinder)
        {
            scene->addChild(builder->createCylinder(info));
            bound.add(info.position);
            info.position += info.dx * 1.5f;
        }

        if (cone)
        {
            scene->addChild(builder->createCone(info));
            bound.add(info.position);
            info.position += info.dx * 1.5f;
        }

        if (capsule)
        {
            scene->addChild(builder->createCapsule(info));
            bound.add(info.position);
        }

        // update the centre and radius to account for all the shapes added so we can position the camera to see them all.
        centre = (bound.min + bound.max) * 0.5;
        radius += vsg::length(bound.max - bound.min) * 0.5;

    }

    if (info.positions)
    {
        FindMatchingBufferInfo findMatches(info.positions);
        scene->accept(findMatches);

        std::cout<<"Found matches "<<findMatches.matches.size()<<" for "<<info.positions<<std::endl;
        for(auto& bufferInfo : findMatches.matches)
        {
            std::cout<<"    "<<bufferInfo<<", "<<bufferInfo->buffer<<", "<<bufferInfo->offset<<", "<<bufferInfo->range<<", "<<bufferInfo->data<<std::endl;
        }
    }

    // create the viewer and assign window(s) to it
    auto viewer = vsg::Viewer::create();

    auto window = vsg::Window::create(windowTraits);
    if (!window)
    {
        std::cout << "Could not create windows." << std::endl;
        return 1;
    }

    viewer->addWindow(window);

    vsg::ref_ptr<vsg::LookAt> lookAt;

    // compute the bounds of the scene graph to help position camera
    //vsg::ComputeBounds computeBounds;
    //scene->accept(computeBounds);
    //vsg::dvec3 centre = (computeBounds.bounds.min + computeBounds.bounds.max) * 0.5;
    //double radius = vsg::length(computeBounds.bounds.max - computeBounds.bounds.min) * 0.6 * 10.0;

    // set up the camera
    lookAt = vsg::LookAt::create(centre + vsg::dvec3(0.0, -radius * 3.5, 0.0), centre, vsg::dvec3(0.0, 0.0, 1.0));

    double nearFarRatio = 0.001;
    auto perspective = vsg::Perspective::create(30.0, static_cast<double>(window->extent2D().width) / static_cast<double>(window->extent2D().height), nearFarRatio * radius, radius * 4.5);

    auto camera = vsg::Camera::create(perspective, lookAt, vsg::ViewportState::create(window->extent2D()));

    // set up the compilation support in builder to allow us to interactively create and compile subgraphs from wtihin the IntersectionHandler
    builder->setup(window, camera->viewportState);

    // add close handler to respond the close window button and pressing escape
    viewer->addEventHandler(vsg::CloseHandler::create(viewer));

    viewer->addEventHandler(vsg::Trackball::create(camera));

    auto memoryBufferPools = vsg::MemoryBufferPools::create("Staging_MemoryBufferPool", window->getDevice(), vsg::BufferPreferences{});
    auto copyBufferCmd = vsg::CopyAndReleaseBuffer::create(memoryBufferPools);

    auto renderGraph = vsg::createRenderGraphForView(window, camera, scene);

    auto commandGraph = vsg::CommandGraph::create(window);
    commandGraph->addChild(copyBufferCmd);
    commandGraph->addChild(renderGraph);

    viewer->assignRecordAndSubmitTaskAndPresentation({commandGraph});

    viewer->compile();

    auto startTime = vsg::clock::now();
    double numFramesCompleted = 0.0;

    vsg::ref_ptr<vsg::vec3Array> positions = info.positions;
    std::set<vsg::ref_ptr<vsg::BufferInfo>> bufferInfos;

    if (info.positions)
    {
        FindMatchingBufferInfo findMatches(info.positions);
        scene->accept(findMatches);
        bufferInfos = findMatches.matches;

        std::cout<<"Found matches "<<findMatches.matches.size()<<" for "<<info.positions<<std::endl;
        for(auto& bufferInfo : findMatches.matches)
        {
            std::cout<<"    "<<bufferInfo<<", "<<bufferInfo->buffer<<", "<<bufferInfo->offset<<", "<<bufferInfo->range<<", "<<bufferInfo->data<<std::endl;
        }
    }

    dynamicUpdates = info.positions.valid() && dynamicUpdates;

    // rendering main loop
    while (viewer->advanceToNextFrame() && (numFrames < 0 || (numFrames--) > 0))
    {
        // pass any events into EventHandlers assigned to the Viewer
        viewer->handleEvents();

        viewer->update();

        if (dynamicUpdates)
        {
            for(auto& v : *positions)
            {
                v.z += 0.01*std::sin(0.01*float(viewer->getFrameStamp()->frameCount)) * v.x;
                //std::cout<<"new vertex "<<v<<std::endl;
            }

            if (frameToWait>0 && waitTimeout>0)
            {
                viewer->waitForFences(frameToWait, waitTimeout);
            }

            // pass the poisitions to GPU
            for(auto& bufferInfo : bufferInfos)
            {
                copyBufferCmd->copy(positions, bufferInfo);
            }
        }

        viewer->recordAndSubmit();

        viewer->present();

        numFramesCompleted += 1.0;
    }

    auto duration = std::chrono::duration<double, std::chrono::seconds::period>(vsg::clock::now() - startTime).count();
    if (numFramesCompleted > 0.0)
    {
        std::cout<<"Average frame rate = "<<(numFramesCompleted / duration)<<std::endl;
    }

    return 0;
}
