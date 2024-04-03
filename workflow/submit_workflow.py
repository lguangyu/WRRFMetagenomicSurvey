#!/usr/bin/env python3

import argparse
import io
import json
import os
import re
import subprocess
import sys


def get_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser()
	ap.add_argument("workflow", type = str, nargs = "?", default = "-",
		help = "workflow config json [stdin]")
	ap.add_argument("--verbose", "-v", action = "store_true",
		help = "increase verbosity [no]")

	# parse and refine args
	args = ap.parse_args()
	if args.workflow == "-":
		args.workflow = sys.stdin

	return args


def get_fp(f, *ka, factory = open, **kw) -> io.IOBase:
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first argument of get_fp() must be str or io.IOBase, "
			"got '%s'" % type(f).__name__)
	return ret


class SlurmRoutineMixin(object):
	"""
	interface routines with slurm
	"""
	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		self.submitted = False
		self.slurm_id = None
		return

	def slurm_submit_job(self, script, sbatch_args = None, verbose = False,
			**kw):
		if sbatch_args is None:
			sbatch_args = list()
		cmd = ["sbatch"] + sbatch_args + [script]
		if verbose:
			print("calling: " + str(cmd), file = sys.stderr)
		slurm_msg = subprocess.check_output(cmd).decode()
		self.submitted = True
		self.slurm_id = self._parse_slurm_submit_id(slurm_msg)
		return

	@staticmethod
	def _parse_slurm_submit_id(s: str) -> str:
		return re.search(r"(\d+)$", s).group(1)


class Target(SlurmRoutineMixin):
	def __init__(self, key, script, mark_as_done = False, dependency = None,
			sbatch_args = None, **kw):
		super().__init__(**kw)
		self.key = key
		self.script = script
		self._check_script_exist()
		self.mark_as_done = mark_as_done
		self.dependency = dependency
		self.dependency_target = list()
		self.sbatch_args = sbatch_args
		# other attibutes
		for k, v in kw.items():
			setattr(self, k, v)
		return

	@property
	def dependency(self):
		return self._dep
	@dependency.setter
	def dependency(self, dep):
		if dep is None:
			self._dep = list()
		elif isinstance(dep, str):
			self._dep = [dep]
		elif isinstance(dep, list):
			self._dep = dep
		else:
			raise TypeError("dependency must be None, str or list, got '%s'"\
				% type(dep).__name__)
		return

	def _check_script_exist(self):
		if not os.path.isfile(self.script):
			raise ValueError("expect script path to be an existing file, but "
				"'%s' seems not" % self.script)
		return

	def slurm_submit_job(self, verbose = False):
		if self.submitted or self.mark_as_done:
			return
		# check and submit dependencies
		for dep in self.dependency_target:
			dep.slurm_submit_job(verbose = verbose)
		# concatenate sbatch arguments
		sbatch_args = list()
		# first, dependency, join ids (presumably sumbitted to this point)
		dep_ids = (":").join([dep.slurm_id for dep in self.dependency_target\
			if not dep.mark_as_done])
		if dep_ids:
			sbatch_args.extend(["--dependency", "afterok:" + dep_ids])
		# add other sbatch args
		if self.sbatch_args:
			sbatch_args.extend(self.sbatch_args)
		super().slurm_submit_job(self.script, sbatch_args, verbose = verbose)
		return


class TargetCollection(object):
	def __init__(self):
		self.target_dict = dict()
		return

	def add_target(self, target: Target):
		self.target_dict[target.key] = target
		return

	def get_target(self, key):
		if not key in self.target_dict:
			raise KeyError("target '%s' does not exist" % key)
		return self.target_dict[key]

	@property
	def target_keys(self):
		return self.target_dict.keys()

	@property
	def targets(self):
		return self.target_dict.values()

	@classmethod
	def from_config_file(cls, f):
		new = cls()
		with get_fp(f, "r") as fp:
			cfg = json.load(fp)
		new._add_all_cfg_target(cfg)
		new._fill_dependency_target()
		return new

	def _add_all_cfg_target(self, cfg):
		script_prefix = cfg["script_prefix"]
		for k, v in cfg["target"].items():
			self._add_cfg_target(k, script_prefix, **v)
		return

	def _add_cfg_target(self, key: str, script_prefix: str = None,
			mark_as_done = False, dependency = None, **kw):
		# assign default variables
		if script_prefix is None: script_prefix = ""
		script = script_prefix + key
		self.add_target(Target(key, script, mark_as_done, dependency, **kw))
		return

	def _fill_dependency_target(self):
		# upon creation, dependency of targets are strings.
		# here to fill another list, dependency_target, with respective Target
		# objects for faster traversal
		for key in self.target_keys:
			target = self.get_target(key)
			target.dependency_target = [self.get_target(i)\
				for i in target.dependency]
		return

	def slurm_submit_all(self, verbose = False):
		# this is a depth-first travesal of the tree (or graph, whatever)
		# at this level, simply iterate every single job
		# recursive nature of target.slurm_submit_job() will hand the traversal
		for target in self.targets:
			target.slurm_submit_job(verbose = verbose)
		return


def main():
	args = get_args()
	tc = TargetCollection.from_config_file(args.workflow)
	tc.slurm_submit_all(verbose = args.verbose)
	return


if __name__ == "__main__":
	main()

