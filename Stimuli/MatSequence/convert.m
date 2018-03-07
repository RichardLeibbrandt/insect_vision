function convert(infile, outfile)
    v = VideoReader(infile);
    vidFrame = readFrame(v);
    out = zeros(size(vidFrame,1), size(vidFrame,2), 3);
    mcframe = rgb2gray(vidFrame);
    out(:,:,1) = mcframe;
    counter = 1;
    while hasFrame(v)
        counter = counter + 1;
        vidFrame = readFrame(v);
        mcframe = rgb2gray(vidFrame);
        out(:,:,counter) = mcframe;
    end
    save(outfile, 'out');
end