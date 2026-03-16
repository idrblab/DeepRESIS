import math

import torch
import torch.nn as nn


class ValueEmbedding(nn.Module):
    def __init__(self, d_feature: int, d_model: int) -> None:
        super().__init__()
        self.value_embedding = nn.Sequential(
            nn.Linear(d_feature, d_model),
            nn.ReLU(),
            nn.Linear(d_model, d_model),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.value_embedding(x)


class PositionalEmbedding(nn.Module):
    def __init__(self, d_model: int, dropout: float, max_len: int = 5000) -> None:
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        pos_emb = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * -(math.log(10000.0) / d_model))
        pos_emb[:, 0::2] = torch.sin(position * div_term)
        pos_emb[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pos_emb", pos_emb.unsqueeze(0))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.pos_emb[:, : x.size(1)].requires_grad_(False)
        return self.dropout(x)


class RNAEmbedding(nn.Module):
    def __init__(self, d_feature: int, d_model: int, dropout: float) -> None:
        super().__init__()
        self.val_emb = ValueEmbedding(d_feature, d_model)
        self.pos_emb = PositionalEmbedding(d_model, dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.pos_emb(self.val_emb(x))


class DrugEmbedding(nn.Module):
    def __init__(self, d_feature1: int, d_feature2: int, d_model: int, dropout: float) -> None:
        super().__init__()
        self.dropout = nn.Dropout(p=0.5)
        self.val_emb1 = ValueEmbedding(d_feature1, d_model)
        self.val_emb2 = ValueEmbedding(d_feature2, d_model)
        self.pos_emb = PositionalEmbedding(d_model, dropout)

    def forward(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        x1 = self.val_emb1(x)
        x2 = self.val_emb2(self.dropout(y))
        return self.pos_emb(torch.cat((x1, x2), 1))


class TransformerEncoder(nn.Module):
    def __init__(self, d_model: int, nhead: int, num_layers: int) -> None:
        super().__init__()
        self.encoder_layer = nn.TransformerEncoderLayer(d_model, nhead, batch_first=True)
        self.encoder = nn.TransformerEncoder(self.encoder_layer, num_layers)

    def forward(self, src: torch.Tensor) -> torch.Tensor:
        return self.encoder(src)


class TransformerDecoder(nn.Module):
    def __init__(self, d_model: int, nhead: int, num_layers: int) -> None:
        super().__init__()
        self.decoder_layer = nn.TransformerDecoderLayer(d_model, nhead, batch_first=True)
        self.decoder = nn.TransformerDecoder(self.decoder_layer, num_layers)

    def forward(self, src: torch.Tensor, tgt: torch.Tensor) -> torch.Tensor:
        return self.decoder(src, tgt)


class Transformer(nn.Module):
    def __init__(
        self,
        dim_1: int,
        dim_2: int,
        dim_3: int,
        d_model: int,
        dropout: float,
        nhead: int,
        num_encoder_layers: int,
        num_decoder_layers: int,
    ) -> None:
        super().__init__()
        self.data_embed_1 = RNAEmbedding(dim_1, d_model, dropout)
        self.data_embed_2 = DrugEmbedding(dim_2, dim_3, d_model, dropout)
        self.encoder_1 = TransformerEncoder(d_model, nhead, num_encoder_layers)
        self.encoder_2 = TransformerEncoder(d_model, nhead, num_encoder_layers)
        self.decoder_1 = TransformerDecoder(d_model, nhead, num_decoder_layers)
        self.decoder_2 = TransformerDecoder(d_model, nhead, num_decoder_layers)

    def forward(self, x1: torch.Tensor, x2: torch.Tensor, x3: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        seq1 = self.data_embed_1(x1)
        seq2 = self.data_embed_2(x2, x3)
        src_1 = self.encoder_1(seq1)
        src_2 = self.encoder_2(seq2)
        out_1 = self.decoder_1(src_1, seq2)
        out_2 = self.decoder_1(src_2, seq1)
        out_2 = out_2.view(out_2.size(0), -1).unsqueeze(1)
        return out_1, out_2


class Classifier(nn.Module):
    def __init__(
        self,
        dim_1: int,
        dim_2: int,
        dim_3: int,
        d_model: int,
        dropout: float,
        nhead: int,
        num_encoder_layers: int,
        num_decoder_layers: int,
    ) -> None:
        super().__init__()
        self.transformer = Transformer(
            dim_1,
            dim_2,
            dim_3,
            d_model,
            dropout,
            nhead,
            num_encoder_layers,
            num_decoder_layers,
        )
        self.fc_x_2 = nn.Linear(d_model * 2, d_model)
        self.fc_out = nn.Linear(d_model * 2, 1)

    def forward(self, x1: torch.Tensor, x2: torch.Tensor, x3: torch.Tensor) -> torch.Tensor:
        x_1, x_2 = self.transformer(x1, x2, x3)
        x_2 = self.fc_x_2(x_2)
        x = self.fc_out(torch.cat((x_1, x_2), 2)).squeeze(1)
        return 1 / (1 + 1 / torch.exp(x))
