
from PIL import Image, ImageFilter, ImageStat


class _Enhce(object):

    def enhce(self, factor):
        """
        Returns an enhanced image.
        :param factor: A floating point value controlling the enhancement.
                       Factor 1.0 always returns a copy of the original image,
                       lower factors mean less color (brightness, contrast,
                       etc), and higher values more. There are no restrictions
                       on this value.
        :rtype: :py:class:`~PIL.Image.Image`
        """
        return Image.blend(self.degenerate, self.image, factor)

class Contrast(_Enhce):
    """Adjust image contrast.
    This class can be used to control the contrast of an image, similar
    to the contrast control on a TV set. An enhancement factor of 0.0
    gives a solid grey image. A factor of 1.0 gives the original image.
    """
    def __init__(self, image,mid):
        self.image = image
        #mean = int(ImageStat.Stat(image.convert("L")).mean[0] + 0.5)
        self.degenerate = Image.new("L", image.size, mid).convert(image.mode)

        if 'A' in image.getbands():
            self.degenerate.putalpha(image.split()[-1])

