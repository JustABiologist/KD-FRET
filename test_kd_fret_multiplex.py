import importlib
import shutil
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


try:
    kd = importlib.import_module("kd_fret_multiplex")
except ImportError as exc:  # pragma: no cover - depends on local scientific env
    kd = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None


def require_helper(name):
    if kd is None:
        raise unittest.SkipTest(f"Could not import kd_fret_multiplex: {IMPORT_ERROR}")
    helper = getattr(kd, name, None)
    if helper is None:
        raise unittest.SkipTest(f"kd_fret_multiplex.{name} is not implemented yet")
    return helper


class TiffSequenceValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp_root = Path.cwd() / "_test_tmp_kd_fret"
        shutil.rmtree(self.tmp_root, ignore_errors=True)
        self.tmp_root.mkdir()

    def tearDown(self):
        shutil.rmtree(self.tmp_root, ignore_errors=True)

    def test_validate_tiff_sequence_frame_count_accepts_exact_indexed_frames(self):
        validate = require_helper("validate_tiff_sequence_frame_count")

        stack_dir = self.tmp_root / "exact"
        stack_dir.mkdir()
        for name in ("xy01use0001.tif", "xy01use0002.tif", "bleach_roi_mask.tif"):
            (stack_dir / name).touch()

        validate(stack_dir, 2, "registered xy01")

    def test_validate_tiff_sequence_frame_count_accepts_zero_based_imagej_frames(self):
        validate = require_helper("validate_tiff_sequence_frame_count")

        stack_dir = self.tmp_root / "zero_based"
        stack_dir.mkdir()
        for name in ("xy01use0000.tif", "xy01use0001.tif"):
            (stack_dir / name).touch()

        validate(stack_dir, 2, "registered xy01")

    def test_validate_tiff_sequence_frame_count_rejects_missing_frames(self):
        validate = require_helper("validate_tiff_sequence_frame_count")

        stack_dir = self.tmp_root / "missing"
        stack_dir.mkdir()
        (stack_dir / "xy01use0001.tif").touch()

        with self.assertRaisesRegex(RuntimeError, "registered xy01"):
            validate(stack_dir, 2, "registered xy01")

    def test_validate_raw_sift_registration_outputs_checks_all_output_sequences(self):
        validate_outputs = require_helper("validate_raw_sift_registration_outputs")
        require_helper("validate_tiff_sequence_frame_count")

        root = self.tmp_root / "outputs"
        root.mkdir()
        measurement = SimpleNamespace(name="measurement1")
        config = SimpleNamespace(
            registered_root=root / "01_registered",
            registered_bleached_root=root / "01_registered_bleached",
            registered_unbleached_root=root / "01_registered_unbleached",
            full_image=False,
        )
        registered_dir = config.registered_root / measurement.name / "xy01"
        registered_dir.mkdir(parents=True)
        (registered_dir / "bleach_roi_mask.tif").touch()
        (registered_dir / "bleach_roi_crop.roi").touch()

        with mock.patch.object(kd, "validate_tiff_sequence_frame_count") as validate_sequence:
            validate_outputs(measurement, config, "xy01", 12)

        expected_calls = [
            mock.call(config.registered_root / measurement.name / "xy01", 12, mock.ANY),
            mock.call(config.registered_bleached_root / measurement.name / "xy01", 12, mock.ANY),
            mock.call(config.registered_unbleached_root / measurement.name / "xy01", 12, mock.ANY),
        ]
        validate_sequence.assert_has_calls(expected_calls, any_order=False)
        self.assertEqual(validate_sequence.call_count, 3)


class RawSiftRetryTests(unittest.TestCase):
    def test_register_and_crop_checked_retries_after_validation_failure(self):
        checked = require_helper("register_and_crop_checked")
        require_helper("validate_raw_sift_registration_outputs")

        measurement = SimpleNamespace(name="measurement1")
        config = SimpleNamespace()
        first_failure = ValueError("wrong frame count")

        with (
            mock.patch.object(kd, "register_and_crop", side_effect=[Path("first"), Path("second")]) as register,
            mock.patch.object(
                kd,
                "validate_raw_sift_registration_outputs",
                side_effect=[first_failure, None],
            ) as validate_outputs,
        ):
            result = checked("ij", measurement, config, "xy01", 12, 1)

        self.assertEqual(result, Path("second"))
        self.assertEqual(register.call_count, 2)
        self.assertEqual(validate_outputs.call_count, 2)
        validate_outputs.assert_called_with(measurement, config, "xy01", 12)

    def test_register_and_crop_checked_raises_after_exhausting_retries(self):
        checked = require_helper("register_and_crop_checked")
        require_helper("validate_raw_sift_registration_outputs")

        measurement = SimpleNamespace(name="measurement1")
        config = SimpleNamespace()

        with (
            mock.patch.object(kd, "register_and_crop", return_value=Path("registered")) as register,
            mock.patch.object(
                kd,
                "validate_raw_sift_registration_outputs",
                side_effect=ValueError("wrong frame count"),
            ) as validate_outputs,
        ):
            with self.assertRaisesRegex(RuntimeError, "wrong frame count"):
                checked("ij", measurement, config, "xy01", 12, 1)

        self.assertEqual(register.call_count, 2)
        self.assertEqual(validate_outputs.call_count, 2)


if __name__ == "__main__":
    unittest.main()
