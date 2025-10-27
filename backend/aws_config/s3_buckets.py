from __future__ import annotations
import yaml
import os
from dotenv import load_dotenv
from botocore import UNSIGNED
from botocore.client import Config
from typing import Any, Dict, Iterator, Tuple


class AccountConfig(dict):
    """Validated configuration for a boto3 client, dict-like with attribute access."""

    VALID_KEYS = {
        "aws_access_key_id",
        "aws_secret_access_key",
        "aws_session_token",
        "region_name",
        "api_version",
        "use_ssl",
        "verify",
        "endpoint_url",
        "config",
    }

    def __init__(self, raw_config: Dict[str, Any]):
        if not isinstance(raw_config, dict):
            raise ValueError("Account config must be a dictionary")

        if all(raw_config.get(attr) for attr in ("aws_access_key_id", "aws_secret_access_key")):
            filtered = {k: v for k, v in raw_config.items() if k in self.VALID_KEYS}
        else:
            filtered = {"config": Config(signature_version=UNSIGNED)}
            if "region_name" in raw_config:
                filtered["region_name"] = raw_config["region_name"]

        super().__init__(filtered)

    def __getattr__(self, item: str) -> Any:
        try:
            return self[item]
        except KeyError:
            raise AttributeError(f"{item} not found in AccountConfig")


class BucketConfig(dict):
    """S3 bucket configuration, dict-like with attribute access."""

    VALID_KEYS = {
        "Bucket",
        "Prefix",
        "Delimiter",
        "MaxKeys",
        "EncodingType",
        "ExpectedBucketOwner",
    }

    def __init__(self, raw_config: Dict[str, Any]):
        if not isinstance(raw_config, dict):
            raise ValueError("Bucket config must be a dictionary")

        filtered = {k: v for k, v in raw_config.items() if k in self.VALID_KEYS}
        super().__init__(filtered)

        if "Bucket" not in self:
            raise ValueError("Bucket config must include a 'Bucket' key")

    def __getattr__(self, item: str) -> Any:
        try:
            return self[item]
        except KeyError:
            raise AttributeError(f"{item} not found in BucketConfig")


class Account:
    """An AWS account containing config and buckets. Dict-like for buckets."""

    def __init__(self, name: str, config: AccountConfig, buckets: Dict[str, BucketConfig]):
        self._name = name
        self.config = config
        self._buckets = buckets
        # Expose buckets as attributes
        for bucket_name, bucket_cfg in buckets.items():
            setattr(self, bucket_name, bucket_cfg)

    def __getitem__(self, key: str) -> BucketConfig:
        return self._buckets[key]

    def __iter__(self) -> Iterator[Tuple[str, BucketConfig]]:
        return iter(self._buckets.items())

    def __len__(self) -> int:
        return len(self._buckets)

    def items(self):
        return self._buckets.items()

    def keys(self):
        return self._buckets.keys()

    def values(self):
        return self._buckets.values()


class S3BucketConfig:
    """Loader for account and bucket configurations from YAML. Dict-like for accounts."""

    def __init__(self, yaml_file: str | None = None):
        load_dotenv(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".env"))
        yaml_file = os.environ.get("S3_BUCKETS_YAML", os.path.join(os.path.dirname(os.path.abspath(__file__)), "s3_buckets.yaml")) if yaml_file is None else yaml_file

        if not os.path.exists(yaml_file):
            raise FileNotFoundError(f"YAML file not found: {yaml_file}")

        with open(yaml_file, "r") as f:
            accounts = yaml.load(f, Loader=yaml.FullLoader)

        if not isinstance(accounts, dict):
            raise ValueError("Invalid YAML format: Expected a dictionary of accounts.")

        self._accounts: Dict[str, Account] = {}

        for account_name, account_config in accounts.items():
            if not isinstance(account_config, dict):
                raise ValueError(f"Invalid configuration for account: {account_name}")

            config_obj = AccountConfig(account_config.get("config", {}) or {})
            buckets_dict = account_config.get("buckets", {}) or {}

            if not isinstance(buckets_dict, dict):
                raise ValueError(f"Invalid buckets block for account: {account_name}")
            if not all(isinstance(v, dict) for v in buckets_dict.values()):
                raise ValueError(f"Invalid bucket config under account: {account_name}")

            bucket_objs = {k: BucketConfig(v) for k, v in buckets_dict.items()}

            account_obj = Account(account_name, config_obj, bucket_objs)
            self._accounts[account_name] = account_obj
            setattr(self, account_name, account_obj)

    def __getitem__(self, key: str) -> Account:
        return self._accounts[key]

    def __iter__(self) -> Iterator[Tuple[str, Account]]:
        return iter(self._accounts.items())

    def __len__(self) -> int:
        return len(self._accounts)

    def items(self):
        return self._accounts.items()

    def keys(self):
        return self._accounts.keys()

    def values(self):
        return self._accounts.values()