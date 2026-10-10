"""Tests for DAQ's spec and speck file writing.

The DAQ is built without __init__, which needs noise spec files, and with a constant signal
and no noise, so every bin of every slice has the same power.
"""

from __future__ import annotations

import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from he6_cres_spec_sims.simulation_blocks.DAQ import DAQ

_FREQ_BINS = 8192
_BINS_PER_CHANNEL = 4096
_HEADER_BYTES = 32
_FOOTER_BYTES = 3
_SLICES = 100
_SLICE_BLOCK = 40  # 100 slices are written as chunks of 40, 40 and 20.
_POWER = 25  # |3 + 4j|**2


def _make_daq(files_dir: Path, *, suffix: str, threshold: float = 20.0, n_acquisitions: int = 1) -> DAQ:
    daq = object.__new__(DAQ)
    daq.config = types.SimpleNamespace(
        daq=types.SimpleNamespace(
            freq_bins=_FREQ_BINS,
            spec_suffix=suffix,
            spec_prefix="t",
            roach_inverted_flag=True,
            acq_length=1.0,
            threshold_factor=1.0,
        )
    )
    daq.thresholds = np.full(_FREQ_BINS, threshold)
    daq.bins = [slice(0, _BINS_PER_CHANNEL), slice(_BINS_PER_CHANNEL, _FREQ_BINS)]
    daq.n_acquisitions = n_acquisitions
    daq.n_channels = 2
    daq.slices_in_roach = _SLICES
    daq.slices_in_spec = _SLICES
    daq.slice_block = _SLICE_BLOCK
    daq.get_signal_array = lambda acq, start_slice, stop_slice: np.full(
        (int(stop_slice - start_slice), _FREQ_BINS), 3 + 4j
    )
    daq.get_noise_array = lambda n_slices: np.zeros((n_slices, _FREQ_BINS))
    daq.spec_file_paths = daq.build_file_paths(n_acquisitions, 2, files_dir)
    return daq


def _write(daq: DAQ, *, max_chunks: int | None = None) -> None:
    with daq._open_spec_files() as spec_files:
        daq._write_acquisitions(spec_files, max_chunks)


def _packet_number(record: bytes) -> int:
    return int.from_bytes(record[9:12], "big")


class SpecFileTest(unittest.TestCase):
    def test_each_slice_is_a_header_with_its_packet_number_then_the_bin_powers(self) -> None:
        record_bytes = _HEADER_BYTES + _BINS_PER_CHANNEL
        with TemporaryDirectory() as tmp:
            daq = _make_daq(Path(tmp), suffix="spec")
            _write(daq)
            for path in daq.spec_file_paths[0]:
                data = path.read_bytes()
                self.assertEqual(len(data), _SLICES * record_bytes)
                for i in range(_SLICES):
                    record = data[i * record_bytes : (i + 1) * record_bytes]
                    self.assertEqual(_packet_number(record), i)
                    self.assertEqual(record[:9] + record[12:_HEADER_BYTES], bytes(29))
                    self.assertEqual(set(record[_HEADER_BYTES:]), {_POWER})


class SpeckFileTest(unittest.TestCase):
    def test_each_slice_lists_every_bin_above_the_threshold(self) -> None:
        triplets = bytes(
            np.array([[j // 256, j % 256, _POWER] for j in range(_BINS_PER_CHANNEL)], dtype=np.uint8).reshape(-1)
        )
        record_bytes = _HEADER_BYTES + len(triplets) + _FOOTER_BYTES
        with TemporaryDirectory() as tmp:
            daq = _make_daq(Path(tmp), suffix="speck", threshold=20.0)
            _write(daq)
            for path in daq.spec_file_paths[0]:
                data = path.read_bytes()
                self.assertEqual(len(data), _SLICES * record_bytes)
                for i in range(_SLICES):
                    record = data[i * record_bytes : (i + 1) * record_bytes]
                    self.assertEqual(_packet_number(record), i)
                    self.assertEqual(record[_HEADER_BYTES:-_FOOTER_BYTES], triplets)
                    self.assertEqual(record[-_FOOTER_BYTES:], bytes(_FOOTER_BYTES))

    def test_bins_below_the_threshold_are_omitted(self) -> None:
        record_bytes = _HEADER_BYTES + _FOOTER_BYTES
        with TemporaryDirectory() as tmp:
            daq = _make_daq(Path(tmp), suffix="speck", threshold=30.0)
            _write(daq)
            for path in daq.spec_file_paths[0]:
                self.assertEqual(path.stat().st_size, _SLICES * record_bytes)


class OpenSpecFilesTest(unittest.TestCase):
    def test_files_are_open_inside_the_context_and_closed_after_it(self) -> None:
        with TemporaryDirectory() as tmp:
            daq = _make_daq(Path(tmp), suffix="spec", n_acquisitions=2)
            with daq._open_spec_files() as spec_files:
                self.assertEqual([len(acq) for acq in spec_files], [2, 2])
                self.assertFalse(any(f.closed for acq in spec_files for f in acq))
            self.assertTrue(all(f.closed for acq in spec_files for f in acq))
            self.assertEqual(
                sorted(p.name for p in Path(tmp).iterdir()), ["t_0_0.spec", "t_0_1.spec", "t_1_0.spec", "t_1_1.spec"]
            )

    def test_existing_content_is_truncated(self) -> None:
        with TemporaryDirectory() as tmp:
            daq = _make_daq(Path(tmp), suffix="spec")
            daq.spec_file_paths[0][0].write_bytes(b"stale")
            with daq._open_spec_files():
                pass
            self.assertEqual(daq.spec_file_paths[0][0].stat().st_size, 0)

    def test_max_chunks_stops_writing_after_that_many_chunks(self) -> None:
        record_bytes = _HEADER_BYTES + _BINS_PER_CHANNEL
        with TemporaryDirectory() as tmp:
            daq = _make_daq(Path(tmp), suffix="spec", n_acquisitions=2)
            _write(daq, max_chunks=1)
            self.assertEqual(daq.spec_file_paths[0][0].stat().st_size, _SLICE_BLOCK * record_bytes)
            self.assertEqual(daq.spec_file_paths[1][0].stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
