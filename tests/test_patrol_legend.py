import pytest

from src.patrol import setLegendMarkerSize


class FakeScatterHandle:
    def __init__(self):
        self.sizes = None

    def set_sizes(self, sizes):
        self.sizes = sizes


class FakeLineHandle:
    def __init__(self):
        self.marker_size = None

    def set_markersize(self, marker_size):
        self.marker_size = marker_size


@pytest.mark.parametrize("attribute", ["legend_handles", "legendHandles"])
def test_legend_marker_size_supports_matplotlib_apis(attribute):
    scatter = FakeScatterHandle()
    line = FakeLineHandle()
    legend = type("FakeLegend", (), {})()
    setattr(legend, attribute, [scatter, line])

    setLegendMarkerSize(legend, marker_size=10)

    assert scatter.sizes == [10]
    assert line.marker_size == 10
